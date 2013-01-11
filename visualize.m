function visualize(room, mic, src, head_angles)

% visualize(room, mic, src, head_angles)
%
% Draw a 3d representation of the room setup and the position of the
% source and listener.  Parameters are the same as those used in
% rlrs.m.  The listener will be plotted as an isosceles triangle, with
% the base going between the two ears and the point in the direction
% the listener is facing.  The source is plotted as a star.

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

if ~exist('head_angles', 'var') head_angles = [0 0 0]; end

ear_ear_nose = [0 -.075 0; 0 .075 0; .1 0 0]';
R = euler_matrix(head_angles);
xyz = mic(:)*[1 1 1] + R * ear_ear_nose;

plot3(src(1), src(2), src(3), '*')
patch(xyz(1,:), xyz(2,:), xyz(3,:), 0);
grid on
axis equal

room = [0 0 0; room(:)'];
axis(room(:))
xlabel('x')
ylabel('y')
zlabel('z')
