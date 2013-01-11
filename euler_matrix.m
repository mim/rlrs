function R = euler_matrix(angles)

% R = euler_matrix(angles)
%
% Create a rotation matrix from a set of three Euler angles.  For more
% information, see the wikipedia page for Euler angles, this is the
% zyz system of rotations, except that the second angle has been
% negated: http://en.wikipedia.org/wiki/Euler_angles 
%
% Angles are measured in radians.  The first angle is approximately
% the rotation of the head in the horizontal plane, with positive
% angles being to the left.  The second angle is approximately the
% elevation of the head, with positive angles being upwards.
% Technically, the first rotation is about the Z axis, the second is
% about the rotated y axis, and the third is again about the Z axis.

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

c = cos(angles);
s = sin(angles);

% zyz
R = [c(1)*c(2)*c(3)-s(1)*s(3), -c(2)*c(3)*s(1)-c(1)*s(3), -c(3)*s(2);
     c(3)*s(1)+c(1)*c(2)*s(3),  c(1)*c(3)-c(2)*s(1)*s(3), -s(2)*s(3);
     c(1)*s(2),                -s(1)*s(2),                 c(2)];
