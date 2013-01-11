=======================================================================
Matlab Recti-Linear Room Simulator
----------------------------------

Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
Distributable under the Gnu General Public License version 3 or later

=======================================================================
Introduction
------------
This code will generate binaural impulse responses from a simulation
of the acoustics of a rectilinear room using the image method.  It has
a number of features that improve the realism and speed of the
simulation.  It can generate a pair of 680 ms impulse responses
sampled at 22050 Hz in 2 minutes on a 2.8 GHz Intel Xeon.  It's easy
to run from within scripts to generate a large set of impulse
responses programmatically.

To improve the realism, it applies anechoic head-related transfer
functions to each incoming reflection, allows fractional delays,
includes frequency-dependent absorption due to walls, includes
frequency- and humidity-dependent absorption due to air, and varies
the speed of sound with temperature.  It also  randomly
perturbs sources in proportion to their distance to the listener to
simulate imperfections in the alignment of the walls.

To improve simulation speed, it performs all calculations in the
frequency domain and the complex exponential generation code is
written in C, it only calculates the Fourier transforms of anechoic
HRTFs as it needs them, and then it caches them, and it culls sources
that are beyond the desired impulse response length or are
significantly quieter than the direct path.


=======================================================================
Contents
--------
These are the functions you might want to use, in order of usefulness:

  rlrs.m
    main function, produces one set of impulse responses

  visualize.m
    takes a subset of the arguments to rlrs and draws a figure showing
    the sizes and orientations of the various parts of the setup

  rt60.m
    takes a supset of the arguments to rlrs and computes the
    reverberation time of the room

  around_rirs.m
    wrapper around rlrs, creates many impulse responses in a
    hemisphere around the listener


=======================================================================
Use
---
As a very simple example of its use, the following code generates an
impulse response that is 15000 samples long, which is 680 ms at a
sampling rate of 22050 Hz.  This captures up to 67th order reflections
in each of the dimensions.  Dimensions are measured in meters, so this
is a relatively large room, with an RT60 estimated to be 350 ms.  The
placement of the microphone and speaker have Gaussian noise with a
standard deviation of 2 cm added to them to simulate errors in
placement.

mex fast_expj.c
room = [9 5 3.5];
mic = [4.5 2.5 1.5]+.02*randn(1,3);
src = [5.5 2.5 1.5]+.02*randn(1,3);
absorp = [.12 .12 .12 .12 .6 .4]';
visualize(room, mic, src)
rir = rlrs(room, mic, src, 15000, absorp);

For a more involved example, the function around_rirs() generates a
set of 169 impulse responses sampled from the upper hemisphere around
a listener's head at a distance of 1 meter.


=======================================================================
See also
--------
It is based in part on two other matlab implementations of the image
method.  The calculations of the mirror source positions are based on
Stephen G. McGovern's rir.m function, available from
http://www.2pi.us/rir.html .  The improvements to realism are based on
Douglas R. Campbell et al.'s roomsim.m function, available from
http://media.paisley.ac.uk/~campbell/Roomsim/ and 
http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5184&objectType=File

This package also includes one set of head-related impulse responses
(HRIRs) from the CIPIC dataset.  This is for subject 21, the KEMAR
dummy with large pinnae.  If you would like to download other HRIRs or
read more about the data, please see
http://interface.cipic.ucdavis.edu/CIL_html/CIL_HRTF_database.htm


=======================================================================
License
-------
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
