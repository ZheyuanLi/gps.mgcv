#include <stdlib.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
int IsAscending (int n, double *x) {
  int flag = 1; double *xi, *xn = x + n - 1;
  for (xi = x; xi < xn; xi++) {
    if (xi[1] <= xi[0]) {flag = 0; break;}
  }
  return flag;
}
SEXP C_IsAscending (SEXP x, SEXP n, SEXP xi) {
  int i = asInteger(xi), l = length(x);
  if (i < 1 || i > l) error("'xi' is out of bound!");
  double *subx = REAL(x) + i - 1;
  int N = asInteger(n);
  if (N > l - i + 1) error("n <= length(x) - xi + 1 required!");
  int flag = IsAscending(N, subx);
  return ScalarInteger(flag);
}
void MakeGrid (double *b, int k, int n, double *x, int rmdup) {
  int ni = k - 1;
  double *ptrb = b, *bend = b + ni, b0, b1;
  double step0 = 1.0 / (n - 1), step;
  int n0 = n - rmdup;
  if (rmdup) x[0] = b[0];
  double *ptrx = x + rmdup, *pend = ptrx + (n0 - 1);
  for (; ptrb < bend; ptrb++, pend += n0) {
    b0 = ptrb[0]; b1 = ptrb[1];
    step = (b1 - b0) * step0;
    if (rmdup) b0 += step;
    for (; ptrx < pend; b0 += step, ptrx++) *ptrx = b0;
    if (b0 < b1) *ptrx++ = b0;
    else {
      if (b1 > 0) step = b1;
      else if (b1 < 0) step = -b1;
      else step = 1.0;
      *ptrx++ = b1 - 2.22e-16 * step;
    }
  }
}
SEXP C_MakeGrid (SEXP b, SEXP n, SEXP rmdup) {
  int K = length(b), N = asInteger(n), RMDUP = asInteger(rmdup);
  SEXP x = PROTECT(allocVector(REALSXP, (K - 1) * (N - RMDUP) + RMDUP));
  MakeGrid(REAL(b), K, N, REAL(x), RMDUP);
  UNPROTECT(1);
  return x;
}
void SmallAtA (int n, double alpha, double *A, double *X) {
  double *A0i, *A0j = A, *A0n = A + n * n, *Aki, *Akj, *Anj;
  double *Xjj = X, *Xij, *Xji;
  double c;
  while (A0j < A0n) {
    Anj = A0j + n;
    Akj = A0j; c = 0.0;
    while (Akj < Anj) {
      c += Akj[0] * Akj[0];
      Akj++;
    }
    c *= alpha; *Xjj = c;
    A0i = A0j + n; Xij = Xjj + 1; Xji = Xjj + n;
    while (A0i < A0n) {
      Aki = A0i; Akj = A0j; c = 0.0;
      while (Akj < Anj) {
        c += Aki[0] * Akj[0];
        Aki++; Akj++;
      }
      c *= alpha; *Xij = c; *Xji = c;
      A0i += n; Xij++; Xji += n;
    }
    A0j += n; Xjj += n + 1;
  }
}
void SmallLtA (int n, double *L, double *A, double *X) {
  double *Xij = X, *Xnn = X + n * n;
  double *Lii, *Lki;
  double *Aij = A, *Akj, *Anj = A;
  double c;
  while (Xij < Xnn) {
    Lii = L; Anj += n;
    while (Aij < Anj) {
      Lki = Lii; Akj = Aij; c = 0.0;
      while (Akj < Anj) {
        c += Lki[0] * Akj[0];
        Lki++; Akj++;
      }
      *Xij++ = c; Lii += n + 1; Aij++;
    }
  }
}
SEXP C_SbarBlocks (SEXP xd, SEXP W, SEXP B) {
  int ord = nrows(W);
  int k1 = length(xd) - 1;
  double *s = REAL(xd);
  double *b = s + k1;
  double *L = REAL(W); int blocksize;
  F77_CALL(dpotf2)("l", &ord, L, &ord, &blocksize FCONE);
  double *Bj = REAL(B);
  blocksize = ord * ord;
  double alpha;
  double *X = malloc(blocksize * sizeof(double));
  SEXP S = PROTECT(alloc3DArray(REALSXP, ord, ord, k1));
  double *Sj = REAL(S);
  while (s < b) {
    SmallLtA(ord, L, Bj, X);
    alpha = 0.5 * (s[1] - s[0]);
    SmallAtA(ord, alpha, X, Sj);
    Bj += blocksize; Sj += blocksize; s++;
  }
  free(X);
  UNPROTECT(1);
  return S;
}
static inline void Block2LTB (int n, double *A, double *L) {
  double *Ajj = A, *Aij, *Anj = A, *Ann = A + n * n;
  double *L0j = L, *Lij;
  while (Ajj < Ann) {
    Aij = Ajj; Lij = L0j; Anj += n;
    while (Aij < Anj) *Lij++ += *Aij++;
    Ajj += n + 1; L0j += n;
  }
}
static inline void ZeroVec (int n, double *x) {
  double *xi = x, *xn = x + n;
  if (n > 0) xi[0] = 0;
  for (xi += n & 1; xi < xn; xi += 2) {
    xi[0] = 0; xi[1] = 0;
  }
}
SEXP C_SbarLTB (SEXP S, SEXP LPBTRF) {
  SEXP Dim = getAttrib(S, R_DimSymbol);
  int *dim = INTEGER(Dim);
  int ord = dim[0];
  int k1 = dim[2];
  int n = k1 + ord - 1;
  SEXP LTB = PROTECT(allocMatrix(REALSXP, ord, n));
  double *L = REAL(LTB), *Lj = L; int blocksize = length(LTB);
  ZeroVec(blocksize, L);
  blocksize = ord * ord;
  double *Sj = REAL(S), *Send = Sj + blocksize * k1;
  while (Sj < Send) {
    Block2LTB(ord, Sj, Lj);
    Sj += blocksize; Lj += ord;
  }
  if (asInteger(LPBTRF)) {
    k1 = ord - 1;
    F77_CALL(dpbtrf)("l", &n, &k1, L, &ord, &blocksize FCONE);
  }
  UNPROTECT(1);
  return LTB;
}
double xtAx (int n, double *A, double *x) {
  double c = 0.0, alpha;
  double *xj = x, *xi, *xn = x + n;
  double *Ajj = A, *Aij;
  while (xj < xn) {
    alpha = xj[0];
    c += Ajj[0] * alpha * alpha;
    Aij = Ajj + 1; xi = xj + 1; alpha += alpha;
    while (xi < xn) {
      c += Aij[0] * xi[0] * alpha;
      Aij++; xi++;
    }
    Ajj += n + 1; xj++;
  }
  return c;
}
void Diff (int n, int k, double *x, double *dx) {
  double *xi = x, *yi = x + k, *xn = x + n, *dxi = dx, alpha, tmp;
  if (k == 1) {
    while (yi < xn) {
      tmp = yi[0] - yi[-1];
      *dxi++ = tmp;
      yi++;
    }
  } else {
    alpha = 1.0 / k;
    while (yi < xn) {
      tmp = (*yi++) - (*xi++);
      tmp *= alpha;
      *dxi++ = tmp;
    }
  }
}
SEXP C_Diff (SEXP x, SEXP k, SEXP n, SEXP xi) {
  int i = asInteger(xi), l = length(x);
  if (i < 1 || i > l) error("'xi' is out of bound!");
  double *subx = REAL(x) + i - 1;
  int N = asInteger(n), K = asInteger(k);
  if (N > l - i + 1) error("n <= length(x) - xi + 1 required!");
  if (K <= 0 || K >= N) error("1 <= k <= n - 1 required!");
  SEXP dx = PROTECT(allocVector(REALSXP, N - K));
  Diff(N, K, subx, REAL(dx));
  UNPROTECT(1);
  return dx;
}
void ComputeLD (double *xt, int k, int d, double *ld) {
  int m = d - 1, p = k - d, i;
  double *dx1, *dx2;
  for (i = 1; i <= m; i++) {
    dx1 = ld + (i - 1) * p; dx2 = dx1 + i;
    while (dx1 < dx2) *dx1++ = 0.0;
    Diff(k - 2 * i, d - i, xt + i, dx1);
  }
}
SEXP C_ComputeLD (SEXP xt, SEXP d) {
  int K = length(xt), D = asInteger(d);
  SEXP ld = PROTECT(allocMatrix(REALSXP, K - D, D - 1));
  ComputeLD(REAL(xt), K, D, REAL(ld));
  UNPROTECT(1);
  return ld;
}
void NullVec (double *ld, int p, int m, double *h) {
  double *dx, *hp = h + p, *hi, c; int j, skip;
  skip = (m - 1); ZeroVec(skip, h);
  hi = h + skip; while (hi < hp) *hi++ = 1.0;
  j = m - 1;
  while (j--) {
    dx = ld + j * p + skip; hi = h + skip; c = 0.0;
    while (hi < hp) {
      c += (*dx) * (*hi);
      *hi = c; hi++; dx++;
    }
  }
  hi = h + skip; c = 0.0;
  while (hi < hp) {c += (*hi) * (*hi); hi++;}
  hi = h + skip; c = 1.0 / sqrt(c);
  while (hi < hp) *hi++ *= c;
}
void NullD (double *ld, int p, int m, double *H) {
  int i; double *h = H;
  for (i = 1; i <= m; i++, h += p) NullVec(ld, p, i, h);
}
SEXP C_NullD (SEXP ld, SEXP m) {
  int P = nrows(ld), M = asInteger(m);
  SEXP H = PROTECT(allocMatrix(REALSXP, P, M));
  NullD(REAL(ld), P, M, REAL(H));
  UNPROTECT(1);
  return H;
}
