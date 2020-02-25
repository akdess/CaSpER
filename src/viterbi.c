////# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

void addVectors(double *, double *, double *, unsigned int);
void setVectorToValue_int(int *, int, unsigned int);
void maxVectorInPlace(double *, int *, double *, unsigned int);

SEXP viterbi(SEXP piZ, SEXP A, SEXP py) {
  PROTECT(piZ = AS_NUMERIC(piZ));
  PROTECT(A = AS_NUMERIC(A));
  PROTECT(py = AS_NUMERIC(py));

  double * prior, * transmat, * obslik;
  prior = NUMERIC_POINTER(piZ);
  transmat = NUMERIC_POINTER(A);
  obslik = NUMERIC_POINTER(py);

  int K, T, tmp;
  K = GET_LENGTH(piZ);

  if (INTEGER(GET_DIM(A))[0] != K || INTEGER(GET_DIM(A))[1] != K) {
    error("The transition matrix must be of size %d x %d", K, K);
  }

  obslik = NUMERIC_POINTER(py);
  T = INTEGER(GET_DIM(py))[1];

  if (INTEGER(GET_DIM(py))[0] != K) {
    error("The observed likelihoods must have %d rows", K);
  }

  SEXP delta_data, psi_data, path_data, d_data, loglik_data;
  double * delta, * d, * loglik;
  int * psi, * path;
  PROTECT(delta_data = allocMatrix(REALSXP, K, T));
  PROTECT(d_data = NEW_NUMERIC(K));
  PROTECT(loglik_data = NEW_NUMERIC(1));
  PROTECT(psi_data = allocMatrix(INTSXP, K, T));
  PROTECT(path_data = NEW_INTEGER(T));
  delta = NUMERIC_POINTER(delta_data);
  d = NUMERIC_POINTER(d_data);
  loglik = NUMERIC_POINTER(loglik_data);
  psi = INTEGER_POINTER(psi_data);
  path = INTEGER_POINTER(path_data);

  int t = 0;
  addVectors(delta + t * K, prior, obslik + t * K, K);
  setVectorToValue_int(psi + t * K, 0, K);

  // FORWARD
  for (int t = 1; t < T; ++t) {
    for (int j = 0; j < K; ++j) {
      addVectors(d, delta + (t - 1) * K, transmat + j * K, K);
      maxVectorInPlace(delta + j + t * K, psi + j + t * K, d, K);
      delta[j + t * K] += obslik[j + t * K];
    }
  }

  // BACKWARD
  t = T - 1;
  maxVectorInPlace(d, path + t, delta + t * K, K); // Use d[0] as temp variable
  loglik[0] = d[0];

  for(t=T-2;t>=0;--t) {
    path[t] = psi[path[t+1] + (t+1)*K];
  }

  SEXP changes_data;
  int * changes;
  PROTECT(changes_data = allocMatrix(INTSXP, 4, T));
  changes = INTEGER_POINTER(changes_data);

  int changesCounter = 0;
  changes[changesCounter + 0 * T] = 0;
  changes[changesCounter + 1 * T] = 0; /* overwritten */
  changes[changesCounter + 2 * T] = path[0];
  changes[changesCounter + 3 * T] = 0;
  changesCounter = 1;

  for(t = 1; t < T; ++t) {
    if (path[t] != path[t - 1]) {
        changes[changesCounter + 0 * T] = t;
        changes[(changesCounter - 1) + 1 * T] = t - 1;
        changes[changesCounter + 2 * T] = path[t];
        changes[changesCounter + 3 * T] = 0;
        changesCounter++;
    }
  }
  changes[(changesCounter - 1) + 1 * T] = T - 1;

  // Reformat to segs
  SEXP segs_data;
  PROTECT(segs_data = allocMatrix(REALSXP, changesCounter, 4));
  double * segs;
  segs = NUMERIC_POINTER(segs_data);
  for (int t = 0; t < changesCounter; ++t) {
    segs[t + 0 * changesCounter] = changes[t + 0 * T] + 1;
    segs[t + 1 * changesCounter] = changes[t + 1 * T] + 1;
    segs[t + 2 * changesCounter] = changes[t + 2 * T] + 1;
    segs[t + 3 * changesCounter] = changes[t + 3 * T];
  }

  // Increment 1 to change from C to R style indexing
  for (int t = 0; t < T; ++t) ++path[t];

  SEXP list, list_names;
  char *names[3] = {"path", "loglik", "seg"};
  PROTECT(list_names = allocVector(STRSXP, 3));
  for (int i = 0; i < 3; ++i) {
    SET_STRING_ELT(list_names, i, mkChar(names[i]));
  }
  PROTECT(list = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(list, 0, path_data);
  SET_VECTOR_ELT(list, 1, loglik_data);
  SET_VECTOR_ELT(list, 2, segs_data);
  setAttrib(list, R_NamesSymbol, list_names);
  UNPROTECT(12);
  return list;

}

void addVectors(double * out, double * u, double * v, unsigned int L) {
  for (unsigned int i=0;i<L;++i)
    out[i] = u[i] + v[i];
  return;
}

void setVectorToValue_int(int * A, int value, unsigned int L) {
  for (unsigned int i = 0; i < L; ++i)
    A[i] = value;
  return;
}

void maxVectorInPlace(
    double * out_value, int * out_index, double * A, unsigned int L) {
  double maxvalue = A[0];
  int index = 0;

  for (unsigned int i = 1; i < L; ++i) {
    if (maxvalue < A[i]) {
      index = i;
      maxvalue = A[i];
    }
  }

  *out_value = maxvalue;
  *out_index = index;
  return;
}
