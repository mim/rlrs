function alpha = atmospheric_attenuation(p,t,r,f)
% alpha = atmospheric_attenuation(pressure, temperature,
%                                 relative_humidity, frequency)
%
% Compute the attenuation due to atmospheric absorption with the
% specified atmospheric conditions.  Alpha is measured in dB/meter,
% pressure in atmospheres, temperature in degrees C, relative_humidity
% in percent (0-100), and frequency in Hz.

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

% Based on http://www.csgnetwork.com/atmossndabsorbcalc.html
% Which seems to come from ISO 9613 part 1, according to
% http://resource.npl.co.uk/acoustics/techguides/absorption/

t = t + 273.15; % convert to kelvin
C = 4.6151-6.8346*((273.16/t).^1.261);
h = r*(10.^C)*p;
tr = t/293.15;  % convert to relative air temperature (re 20 deg C)
frO = p*(24+4.04e4*h*(0.02+h)/(0.391+h));
frN = p*(tr.^-0.5)*(9+280*h*exp(-4.17*((tr.^(-1/3))-1)));
alpha = 8.686*f.^2.*(1.84e-11*(1/p)*sqrt(tr) + ...
                     (tr.^-2.5)*(0.01275*(exp(-2239.1/t)./(frO+f.^2/frO)) + ...
                                 0.1068*(exp(-3352/t)./(frN+f.^2/frN))));

% Avoid giving back 0
alpha = alpha + eps;
