function ea = euler_angles_for_look(vDst)

% ea = euler_angles_for_look(vDst)
% 
% Calculate the Euler rotation angles to go from looking down the x-axis to
% looking in direction vDst. The rotation will only involve azimuth and
% elevation, the third angle will always be 0.  This will only work with
% the coordinates that the euler_matrix function requires and will only
% work for rotating the x-axis.
%
% To test this:
% assert(all(vDst == [1 0 0] * euler_matrix(euler_angles_for_look(vDst))))

% Copyright (C) 2013 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

vSrc = [1;0;0];
vDst = vDst(:);

ea = [0 0 0];
a = angle(vSrc(1:2), vDst(1:2));
ea(1) = -a;
P = [cos(a) sin(a) 0; 0 0 1];
%P = [-sin(ea(1)) cos(ea(1)) 0; 0 0 1];
ea(2) = -angle(P * vSrc, P * vDst);


function a = angle(v1, v2)
mag = acos(v1' * v2 ./ (norm(v1) * norm(v2)));
sgn = v1(1)*v2(2) - v1(2)*v2(1);
a = sign(sgn) * mag;
