////# Author: Daniel Lai https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "hmmcopy.h"
         
//Register "forward_backward" and "viterbi"
static const R_CallMethodDef callMethods[]  = {
	{"forward_backward", (DL_FUNC) &forward_backward, 3},
	{"viterbi", (DL_FUNC) &viterbi, 3},
	{NULL, NULL, 0}
};

void R_init_HMMcopy(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
