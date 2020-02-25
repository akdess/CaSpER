//# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

double normalizeInPlace(double *, unsigned int);
void multiplyInPlace(double *, double *, double *, unsigned int);
void multiplyMatrixInPlace(double *, double *, double *, unsigned int);
void transposeSquareInPlace(double *, double *, unsigned int);
void outerProductUVInPlace(double *, double *, double *, unsigned int);
void componentVectorMultiplyInPlace(double *, double *, double *, unsigned int);

SEXP forward_backward(SEXP piZ, SEXP A, SEXP py) {
  // Smarter memory management definitely possible
  // Consider freeing objects sooner
  // Allocate objects in order of decreasing size

  PROTECT(piZ = AS_NUMERIC(piZ));  // K
  PROTECT(A = AS_NUMERIC(A));      // K, K
  PROTECT(py = AS_NUMERIC(py));    // K, T

  double *init_state_distrib, *transmat, *obslik, *transmatT;
  init_state_distrib = NUMERIC_POINTER(piZ);
  transmat = NUMERIC_POINTER(A);
  obslik = NUMERIC_POINTER(py);

  int K, T, tmp;

  K = GET_LENGTH(piZ);

  SEXP transmatT_data;
  PROTECT(transmatT_data = NEW_NUMERIC(K * K));
  transmatT = NUMERIC_POINTER(transmatT_data);
  transposeSquareInPlace(transmatT, transmat, K);

  if (INTEGER(GET_DIM(A))[0] != K || INTEGER(GET_DIM(A))[1] != K) {
    error("The transition matrix must be of size %d x %d", K, K);
  }

  obslik = NUMERIC_POINTER(py);
  T = INTEGER(GET_DIM(py))[1];

  if (INTEGER(GET_DIM(py))[0] != K) {
    error("The observed likelihoods must have %d rows", K);
  }

  SEXP scale_data, alpha_data, beta_data, gamma_data, loglik_data;
  double * scale, * alpha, * beta, * gamma, * loglik;
  PROTECT(scale_data = NEW_NUMERIC(T));
  PROTECT(alpha_data = allocMatrix(REALSXP, K, T));
  PROTECT(beta_data = allocMatrix(REALSXP, K, T));
  PROTECT(gamma_data = allocMatrix(REALSXP, K, T));
  PROTECT(loglik_data = NEW_NUMERIC(1));
  scale = NUMERIC_POINTER(scale_data);
  alpha = NUMERIC_POINTER(alpha_data);
  beta = NUMERIC_POINTER(beta_data);
  gamma = NUMERIC_POINTER(gamma_data);
  loglik = NUMERIC_POINTER(loglik_data);

  // *********FOWARD**********
  int t = 0;
  multiplyInPlace(alpha + t * K, init_state_distrib, obslik + t * K, K);
  scale[t] = normalizeInPlace(alpha + t * K, K);

  SEXP m_data;
  double *m;
  PROTECT(m_data = NEW_NUMERIC(K));
  m = NUMERIC_POINTER(m_data);

  for (t = 1; t < T; ++t) {
    multiplyMatrixInPlace(m, transmatT, alpha + (t - 1) * K, K);
    multiplyInPlace(alpha + t * K, m, obslik + t * K, K);
    scale[t] = normalizeInPlace(alpha + t * K, K);
  }

  loglik[0] = 0;
  for (t = 0; t < T; ++t) {
    loglik[0] += log(scale[t]);
  }

  // *********BACKWARD**********

  t = T - 1;
  for(int d = 0; d < K; ++d) {
    beta[d + t * K] = 1;
    gamma[d + t * K] = alpha[d + t * K];
  }

  double *eta, *b, *squareSpace;

  SEXP b_data, eta_data, squareSpace_data, dim;
  PROTECT(b_data = NEW_NUMERIC(K));
  PROTECT(eta_data = NEW_NUMERIC(K * K * T));
  PROTECT(squareSpace_data = allocMatrix(REALSXP, K, K));
  b = NUMERIC_POINTER(b_data);
  eta = NUMERIC_POINTER(eta_data);
  squareSpace = NUMERIC_POINTER(squareSpace_data);

  PROTECT(dim = NEW_INTEGER(3));
  INTEGER(dim)[0] = K; INTEGER(dim)[1] = K; INTEGER(dim)[2] = T;
  SET_DIM(eta_data, dim);

  for(int d = 0; d < (K * K); ++d) {
    *(eta + d + t * K * K) = 0;
  }

  // We have to remember that the 1:T range in R is 0:(T-1) in C
  for(t = (T - 2); t >= 0; --t) {

    // setting beta
    multiplyInPlace(b, beta + (t + 1) * K, obslik + (t + 1) * K, K);
    multiplyMatrixInPlace(m, transmat, b, K);
    normalizeInPlace(m, K);
    for(int d = 0; d < K; ++d) {
      beta[d + t * K] = m[d];
    }
    // using "m" again as valueholder

    // setting eta, whether we want it or not in the output
    outerProductUVInPlace(squareSpace, alpha + t * K, b, K);
    componentVectorMultiplyInPlace(
      eta + t * K * K, transmat, squareSpace, K * K);
    normalizeInPlace(eta + t * K * K, K * K);

    // Setting gamma
    multiplyInPlace(m, alpha + t * K, beta + t * K, K);
    normalizeInPlace(m, K);
    for(int d = 0; d < K; ++d) {
      gamma[d + t * K] = m[d];
    }
  }

  SEXP list, list_names;
  char *names[5] = {"rho", "alpha", "beta", "xi", "loglik"};
  PROTECT(list_names = allocVector(STRSXP, 5));
  for (int i = 0; i < 5; ++i) {
    SET_STRING_ELT(list_names, i, mkChar(names[i]));
  }
  PROTECT(list = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(list, 0, gamma_data);
  SET_VECTOR_ELT(list, 1, alpha_data);
  SET_VECTOR_ELT(list, 2, beta_data);
  SET_VECTOR_ELT(list, 3, eta_data);
  SET_VECTOR_ELT(list, 4, loglik_data);
  setAttrib(list, R_NamesSymbol, list_names);
  UNPROTECT(16);
  return list;
}

// returns the normalization constant used.
double normalizeInPlace(double * A, unsigned int N) {
  double sum = 0;

  for(unsigned int n = 0; n < N; ++n) {
    sum += A[n];
    if (A[n] < 0) {
      error("Cannot normalize a vector with negative values.");
    }
  }

  if (sum == 0) {
    error("Will not normalize a vector of all zeros");
  } else {
    for(unsigned int n = 0; n < N; ++n) {
      A[n] /= sum;
    }
  }
  return sum;
}

void multiplyInPlace(double * result, double * u, double * v, unsigned int K) {
  for(unsigned int n = 0; n < K; ++n)
    result[n] = u[n] * v[n];
  return;
}

void multiplyMatrixInPlace(
    double * result, double * trans, double * v, unsigned int K) {
  for(unsigned int d = 0; d < K; ++d) {
    result[d] = 0;
    for (unsigned int i = 0; i < K; ++i){
      result[d] += trans[d + i * K] * v[i];
    }
  }
  return;
}

void transposeSquareInPlace(double * out, double * in, unsigned int K) {
  for(unsigned int i = 0; i < K; ++i){
    for(unsigned int j = 0; j < K; ++j){
      out[j + i * K] = in[i + j * K];
    }
  }
  return;
}

void outerProductUVInPlace(
    double * out, double * u, double * v, unsigned int K) {
  for(unsigned int i = 0; i < K; ++i){
    for(unsigned int j = 0; j < K; ++j){
      out[i + j * K] = u[i] * v[j];
    }
  }
  return;
}

// this works for matrices also if you just set the length "L" to be the
// right value, often K*K, instead of just K in the case of vectors
void componentVectorMultiplyInPlace(
    double * out, double * u, double * v, unsigned int L) {
  for(unsigned int i = 0; i < L; ++i)
	  out[i] = u[i] * v[i];
  return;
}
