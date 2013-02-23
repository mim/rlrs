function around_rirs(savefile, varargin)

% around_rirs(savefile, [key1, value1, ...])
%
% Simulate impulse responses for sources coming from different azimuths and
% elevations around a stationary listener.  Sources are jiggled around the
% position they are supposed to be in with gaussian noise.  The look
% direction (option 'look_vec') determines the coordinates and axes that
% azimuth and elevation are computed in.  So azimuth of 0 and elevation of
% 0 corresponds to a source directly in the look direction at a distance of
% 'dist'.
%
%
% Options:
% azs       ([-90:15:90]) azimuths to compute IRs for
% els       ([0:15:180]) elevations to compute IRs for
% look_vec  ([1 0 0]) direction listener is looking, az and el are relative to this
% reps      (1) number of repetitions per (az,el)
% pos_std   (.02) stddev of noise in positioning mic and srcs per rep
% room      ([9 5 3.5]) size of the room in meters
% head_loc  ([4.5 2.5 1.5]) location of head before noise is added
% dist      (1) distance from mic to src in meters
% absorp    ([.1 .1 .1 .6 .4]') absorption of room walls (see rlrs)
% ir_length (.68) length of the computed impulse response (sec)
% sr        (22050) sampling rate of computed impulse responses
% jitter    (1e-3) fraction of their distance to jiggle virtual srcs

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

[pos_std, azs, els, dist, reps, ir_length, sr, ...
 jitter, absorp, room, center, look_dir] = process_options(varargin, ...
        'pos_std', .02, 'azs', [-90:15:90], 'els', [0:15:180], ...
        'dist', 1, 'reps', 1, 'ir_length', 15001/22050, 'sr', 22050, ...
        'jitter', 1e-3, 'absorp', [.12 .12 .12 .12 .6 .4]', ...
        'room', [9 5 3.5], 'head_loc', [4.5 2.5 1.5], 'look_dir', [1 0 0])

% Convert to radians
azs = azs * pi/180;
els = els * pi/180;

% Convert to samples
ir_length = round(ir_length * sr / 4) * 4

mic = center + pos_std*randn(size(center));
rotation = euler_matrix_for_look(look_dir);

for eli=1:length(els)
  for azi=1:length(azs)
    % Convert azimuth, elevation, distance into xyz, add some
    % positioning error
    el = els(eli);
    az = azs(azi);
    src_offset = dist*[cos(az)*cos(el) sin(az) cos(az)*sin(el)];
    src_offset = src_offset * rotation;

    for rep=1:reps
      src = center + src_offset + pos_std*randn(size(center));
      srcs{azi,eli,rep} = src;
      
      brir{azi,eli,rep} = rlrs(room, mic, src, ir_length, absorp, ...
                               'sr', sr, 'pos_std', jitter, ...
                               'look_dir', look_dir);
    end
    
    save(savefile);
  end  
end
