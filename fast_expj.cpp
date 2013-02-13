/*
 * y = fast_expj(seed, N);
 *
 * Generate a conjugate symmetric complex sinusoid from a particular
 * seed element, i.e. y(n) = seed^n where seed = exp(-j w0).  The
 * first argument is the seed, the second is the length of the output,
 * both must be doubles.  The symmetry is the same as when you use
 * matlab's fft function, i.e. n = [0 1 2 ... N/2 -N/2 -N/2+1 ... -2
 * -1].
 *
 * Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
 * Distributable under the GPL version 3 or higher
 */
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{ 
  if(nrhs != 2)
    mexErrMsgTxt("Usage: y = fast_expj(seed, N)");
  if(!(mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])))
    mexErrMsgTxt("Arguments must be doubles");
  if(!mxIsComplex(prhs[0]))
    mexErrMsgTxt("Seed must be complex");

  double seed_re = mxGetPr(prhs[0])[0];
  double seed_im = mxGetPi(prhs[0])[0];
  int N = mxGetPr(prhs[1])[0];
  int mid = N/2;

  /* plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
  double *re = mxGetPr(plhs[0]);
  double *im = mxGetPi(plhs[0]); */

  plhs[0] = mxCreateNumericMatrix(N, 1, mxSINGLE_CLASS, mxCOMPLEX);
  float *re = (float *)mxGetPr(plhs[0]);
  float *im = (float *)mxGetPi(plhs[0]);
  
  re[0] = 1;
  im[0] = 0;
  
  /* When the delay is integral, the copying is identical to looping
     all the way to i = N-1.  When the delay is not integral, however,
     looping all the way to i=N-1 does something like the hilbert
     transform to the higher frequencies, giving them all a
     frequency-independent phase offset.  Copying the conjugate
     symmetric version from across the origin avoids this problem and
     is the right thing to do. */
  int i;
/*   for (i = 1; i <= N; i++) { */
  for (i = 1; i <= mid; i++) {
    re[i] = re[i-1]*seed_re - im[i-1]*seed_im;
    im[i] = re[i-1]*seed_im + im[i-1]*seed_re;

    /* when N is even, mid = N-mid, so this will assign this sample to
       the negative frequencies. doesn't affect the real part of the
       fourier transform at all.*/
    re[N-i] = re[i];
    im[N-i] = -im[i];
  }

  /* matlab prefers it when this is flipped back */
  if((N % 2) == 0)
    im[mid] = -im[mid]; 
}
