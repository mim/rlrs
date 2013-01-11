function [imps,rooms] = random_rooms(varargin)

% [imps,rooms] = random_rooms([key1, value1, ...])
%
% Make impulse responses from a number of random positions in random
% rooms.  Options provide bounds on randomness, generally everything
% is uniformly distributed.  First a random rectilinear room is
% created (R_room) with a random absorption characteristic (R_absorp).
% Then the listener is placed in the room randomly at least a certain
% distance from the walls (R_mic).  Then the source is placed
% randomly, but close to the listener, so that it is generally within
% 1-2 meters and a certain distance from the walls (R_src).  The
% listener faces the source within the tolerance of R_angle.
%
% Options are:
% N_imp        (1000): number of rooms
% ir_length    (1024): length of impulse responses in samples
% R_room (2x3 matrix): bounds on room dimensions
% R_absorp  (2x6 mat): bounds on absorption coefficients
% R_angle   (2x1 mat): bounds on angle deviation from src
% R_mic     (2x3 mat): bounds on mic placement, relative to room
% R_src     (2x3 mat): bounds on source placement, relative to room
% seed            (0): seed for random number generators
% rlrs_opts      ({}): extra options to pass to rlrs directly

N_imp = 1000;
ir_length = 1024;
R_room = [4 4 2.5; 10 10 6];
%R_absorp = [.1 .1 .1 .1 .5 .3; .15 .15 .15 .15 .7 .5];
R_absorp = [.1 .1 .1 .1 .1 .2; .15 .15 .15 .15 .15 .4];
R_angle = [-5*pi/180; 5*pi/180];

% x=0, y=0, z=0; x=Lx, y=Ly, z=Lz
R_mic = [1 1 1; -1 -1 2];
R_src = [1 1 1; -1 -1 2];

seed = 0;

% options that get passed straight through to rlrs
rlrs_opts = {};

% Change default options if the caller wants
[N_imp, ir_length, R_room, R_absorp, R_mic, R_src, seed, R_angle, ...
 rlrs_opts] = process_options(varargin, 'N_imp', N_imp, 'ir_length', ...
                              ir_length, 'R_room', R_room, 'R_absorp', ...
                              R_absorp, 'R_mic', R_mic, 'R_src', ...
                              R_src, 'seed', seed, 'R_angle', R_angle, ...
                              'rlrs_opts', rlrs_opts);

% Make reproducible
rand('seed', seed); randn('seed', seed);

% Draw all room parameters first, so everything is reproducible
% regardless of randomness in the actual impulse response
% generation
disp('Drawing parameters...')
for i=1:N_imp
  r.room = draw_params(R_room);
  r.absorp = draw_params(R_absorp);
  r.mic = draw_params(R_mic, r.room);

  r.src = [-1 -1 -1];
  [scale, offset] = range2so(R_src, r.room);
  attempt=1;
  while any(r.src < offset) || any(r.src > offset+scale)
    if (attempt < 10)
      % Put the sources close to the mics, generally about 2 meters away
      xy_angle = rand(1)*2*pi;
      radius = random('gamma', 2, .5) + 1;
      r.src = r.mic + [radius*cos(xy_angle) radius*sin(xy_angle) ...
                       .1*rand(1)];
    else
      r.src = draw_params(R_src, r.room);
    end
    attempt = attempt + 1;
  end
  % Have the listener looking approximately at the source
  xy_angle = atan2(r.src(2)-r.mic(2), r.src(1)-r.mic(1));
  r.angle = xy_angle + R_angle(1) + rand(1)*diff(R_angle);

  rooms(i) = r;
end


% Visualize the mic-source relationships
S = cat(1,rooms.src); M = cat(1, rooms.mic);
plot3([S(:,1) M(:,1)]', [S(:,2) M(:,2)]', [S(:,3) M(:,3)]', ...
      M(:,1), M(:,2), M(:,3), '.')
axis equal
axis([0 R_room(2,1) 0 R_room(2,2) 0 R_room(2,3)]);
grid on, drawnow

% Run the room simulations
disp('Simulating rooms...')
imps = zeros(ir_length, 2, N_imp);
for i=1:N_imp
  r = rooms(i);
  imps(:,:,i) = rlrs(r.room, r.mic, r.src, ir_length, r.absorp, ...
                     'head_angles', [-r.angle 0 0], rlrs_opts{:});
end



function p = draw_params(varargin)
% Draw random parameters in the given range.  

[scale, offset] = range2so(varargin{:});
p = rand(size(scale)).*scale + offset;


function [scale, offset] = range2so(range, upper)
% Compute the scale and offset of a hype-rectangle from a range and
% upper bound.  Range is a 2xN matrix, specifying the min and max of
% each of N dimensions to be drawn.  Upper is 1xN.  Positive entries
% in range are left as is, negative entries are added to upper.  Lower
% bounds are added to 0, upper bounds are subtracted from the values
% in upper.

neg = range < 0;

if nargin < 2
  if any(neg(:))
    error('For negative range, you must supply an upper limit')
  else
    upper = 0;
  end
end

lims = (1-neg).*range + neg.*([range(1,:); upper+range(2,:)]);
offset = lims(1,:);
scale = diff(lims,1);
