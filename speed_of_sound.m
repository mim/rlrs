function c = speed_of_sound(temperature)

% Calculate the speed of sound in air of a specified temperature.
% The speed is measured in meters per second, the temperature in
% degrees celsius.

% Copyright (C) 2003  Douglas R Campbell
% Distributable under the GPL version 3 or higher

c = round(331*sqrt(1+0.0036*temperature));
