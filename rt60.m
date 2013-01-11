function T = rt60(room,absorp,temperature,humidity)

% T = rt60(room,absorp,temperature,humidity)
% 
% Norris-Eyring estimate of the reverberation time of a room with
% dimensions specified in the 3-vector room.  Absorp is the same as
% the one passed to rlrs, either 1x1, 1x6, 6x1, or 6x6.  Rows specify
% the wall, columns the frequency.  Reverberation time is the RT60,
% the time it takes for a sound to decay by 60 dB, measured in
% seconds.

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

if nargin < 4, humidity = 50; end
if nargin < 3, temperature = 20; end
if size(absorp, 1) == 1, absorp = repmat(absorp, 6, 1); end

c = speed_of_sound(temperature);
volume = prod(room);
wall_areas = room([2 2 3 3 1 1]) .* room([3 3 1 1 2 2]);
surface_area = sum(wall_areas);

effective_absorbing_area = wall_areas * absorp;
mean_wall_absorption = effective_absorbing_area / surface_area;

frequencies = [125 250 500 1000 2000 4000];
air_absorption = atmospheric_attenuation(1, temperature, ...
                                         humidity, frequencies);
if size(absorp,2) == 1
  air_absorption = mean(air_absorption);
end

% % Sabine estimate
% T = (55.25/c) * volume ./ effective_absorbing_area;

% Norris Eyring
T = (55.25/c) * volume ./ ...
    (4*air_absorption*volume - ...
     surface_area*log(1 - mean_wall_absorption+eps));
