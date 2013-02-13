function y = fast_expj(seed, N)

% y = fast_expj(seed, N);
% 
% Generate a conjugate symmetric complex sinusoid from a particular
% seed element, i.e. y(n) = seed^n where seed = exp(-j w0).  The first
% argument is the seed, the second is the length of the output, both
% must be doubles.  The symmetry is the same as when you use matlab's
% fft function, i.e. n = [0 1 2 ... N/2 -N/2 -N/2+1 ... -2 -1].

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

error(['Please compile the C version of this function by running ' ...
       '"mex fast_expj.cpp"'])
