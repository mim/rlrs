function rir = rlrs(room, mic, source, N, absorp, varargin)
% rir = rlrs(room, mic, source, samples, absorp, [key1, value1, ...])
%
% Simulate the acoustics of a rectilinear room.  Has many of the
% features of Roomsim, with the simplicity of rir.m.  Uses
% direction-dependent CIPIC HRIRs for each echo, frequency-dependent
% wall responses (using table from Roomsim), and temperature-dependent
% speed of sound.  Performs calculations in the frequency domain.
%
% ROOM is a 3-vector describing the x, y, and z dimensions of the room
% in meters.  MIC and SOURCE are the xyz coordinates of the mic and
% source, in meters.  SAMPLES is the length of the impulse response in
% samples.
%
% ABSORP is a matrix of absorption coefficients of the walls and can
% be 1x1, 1x6, 6x1, or 6x6.  Each row is one of the walls and each
% column is a frequency.  Walls are x=0, x=Lx, y=0, y=Ly, z=0, z=Lz.
% Frequencies are [125, 250, 500, 1000, 2000, 4000] Hz.  If only one
% row is specified, then all of the walls get that response.  If only
% one column is specified, then absorption is frequency independent
% (and the calculation is much faster).
%
% Other arguments that can be passed in as 'key', value pairs:
% 'sr' is the sampling rate, defaults to 22050
% 'hrir_file' is the file to read CIPIC-style head-related impulse
%    responses from, it defaults to 'subject_021_hrir'
% 'head_angles', is a set of three Euler angles specifying the rotation
%    of the head relative to looking down the positive x-axis.  The three
%    angles are basically azimuth, elevation, and tilt, but see
%    euler_matrix.m for specifics.  
% 'thresh' is the lower limit of amplitude for sources to be
%    included, in dB relative to the direct path, defaults to 120
% 'humidity' is the relative humidity of the air, defaults to 65
% 'temperature' is in degrees C
% 'pos_std' is the fraction of the distance to the mics to perturb
%    virtual sources, defaults to 1e-3.  
% 'frac_delay' is a flag indicating whether fractional delays
%   should be allowed, defaults to 0

% Copyright (C) 2008 Michael Mandel <mim at ee columbia edu>
% Distributable under the GPL version 3 or higher

%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%
[humidity, temperature, thresh, pos_std, hrir_file, head_angles, sr, ...
 frac_delay, verbose] = ...
    process_options( varargin, 'humidity', 65, 'temperature', ...
                     20, 'thresh', 120, 'pos_std', 1e-3, ...
                     'hrir_file', 'subject_021_hrir', ...
                     'head_angles', [0 0 0]', 'sr', 22050, ...
                     'frac_delay', 0, 'verbose', 0);

tic

c = speed_of_sound(temperature);
N = ceil(N/2)*2;  % Make N even

% Estimate maximum number of reflections, includes lots of extra
% virtual sources in other dimensions that will be culled out later
n = ceil(N/min(room) * c/sr);
if verbose, fprintf('order=%d, ', n); end

[head_width, els, azs, hrtf_cache, N, time_margin, pos3d, inds] = ...
    load_hrtf_data(hrir_file, N, sr);
% Note that N is now incremented by time_margin

[absorp, air_absorp] = ...
    setup_absorp(absorp, N, sr, temperature, humidity);

%%%%%%%%%%%%%%%%%% Position Calculation %%%%%%%%%%%%%%%%%%%%
% Calculate positions of virtual sources and distances to them
% From rir: http://www.2pi.us/rir.html which includes equations
nn   = single(-n:1:n);                  % Index for the sequence
rms  = nn+0.5-0.5*(-1).^nn;             % Part of equations 2,3,& 4
srcs = (-1).^(nn);                      % part of equations 2,3,& 4
xi   = srcs*source(1)+rms*room(1)-mic(1);    % Equation 2 
yj   = srcs*source(2)+rms*room(2)-mic(2);    % Equation 3 
zk   = srcs*source(3)+rms*room(3)-mic(3);    % Equation 4 

[x,y,z] = meshgrid(xi,yj,zk);           % convert vectors to 3D matrices
x = x(:); y = y(:); z = z(:);

[d, time, x, y, z] = keep_in_range(x, y, z, sr, c, N, time_margin, ...
                                   frac_delay, verbose);

% Randomize positions of farther-away sources
x = x + pos_std*d .* randn(size(d));
y = y + pos_std*d .* randn(size(d));
z = z + pos_std*d .* randn(size(d));
[d, time, x, y, z] = keep_in_range(x, y, z, sr, c, N, time_margin, ...
                                   frac_delay, verbose);


%%%%%%%%%%%%%%%%%% Derived Quantities %%%%%%%%%%%%%%%%%%%%%%
[azi,eli] = closest_2d_index(pos3d, head_angles, [x y z], inds);

% Calculate a 6xN matrix of reflection counts off each wall
refl = [room_num_to_refl_count(floor(x' / room(1)));
        room_num_to_refl_count(floor(y' / room(2)));
        room_num_to_refl_count(floor(z' / room(3)))];

log_refl = log(1-absorp+eps)';
log_d    = single(log(d));
nat2dB   = single(20 / log(10));
min_logd = -nat2dB*min(log_d);
seeds    = double(exp(-j*2*pi*time/N));

%%%%%%%%%%%%%%%%%%%%% Combination %%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through sources, adding in effects of everything
FR  = single(zeros(N,2));
kept = zeros(size(d));
test_point = ceil(size(log_refl, 1) / 4);
for vs=1:length(d)
  % Wall absorption and 1/r shrinkage.
  log_ampl = log_refl * refl(:,vs) - log_d(vs);
  % Air absorption, esp high freqs
  log_ampl = log_ampl - d(vs)*air_absorp;

  % Ignore sources more than thresh dB below the direct path to
  % save computation
  if nat2dB*log_ampl(test_point) < min_logd - thresh
    % TODO: keep track of these and count them in aggregate
    % sum up 1/refl(:,vs) for each delay or set of delays
    continue
  end
  kept(vs) = 1;

  % Use a mex function to generate the complex sinusoid
  % iteratively as exp(j w0) .^ n
  src = fast_expj(seeds(vs), N);
  src = src .* exp(log_ampl);
  
  % HRTF filtering
  [hrtf, hrtf_cache] = lazy_hrtf(azi(vs), eli(vs), hrtf_cache);
  src = repmat(src,1,2) .* hrtf;
  FR = FR + src;
end

if verbose
  fprintf('kept %d sources = %f per sample\n', sum(kept), sum(kept) / N);
end
% k = find(kept);
% plot3(x(k), y(k), z(k), '.'), grid on, axis equal

rir = double(real(ifft(FR)));
rir = rir(1:end-time_margin,:);

if frac_delay && 0
  % Need to clean up high frequency ringing issues
  [B,A] = butter(2, 0.99);
  rir = filter(B, A, rir);
end

if verbose, toc, end


%%%%%%%%%%%%%%%%%%%%%%%% Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d, time, x, y, z] = ...
    keep_in_range(x, y, z, sr, c, N, time_margin, frac_delay, verbose)
% Calculate distances, azimuths, and elevations in cipic terms
d = sqrt(x.^2 + y.^2 + z.^2);
if frac_delay
  time = sr*d/c + 1;
else
  time = round(sr*d/c) + 1;
end

% Cull sources that are too far away
keep    = time < N - time_margin;
if verbose, 
  fprintf('keeping %d of %d reflections\n', sum(keep), length(keep))
end

d = d(keep); time = time(keep);
x = x(keep); y = y(keep); z = z(keep);


function cnt = room_num_to_refl_count(room_num)
% From a vector of virtual room number indices (from any dimension),
% compute a 2xN matrix of reflection counts.  The first row is the
% number of reflections at the wall that pass through the origin, the
% second is the number of reflections that pass through the other wall
% parallel to it.  Equations 9.X in rir.html.

cnt = [abs(.5*room_num -.25 + .25*(-1).^room_num);
       abs(.5*room_num +.25 - .25*(-1).^room_num)];


function ind = closest_index(target, x)
% Find the index in target of the closest element to each element
% in x.  Both target and x are vectors.

d = dist(target(:), x(:));
ind = argmin(d, 1);


function [azi,eli] = closest_2d_index(target, head_angles, x, inds)
% Find the indices of the azimuth and elevation of the closest target
% point in 3-space to each point in x in 3-space.  Target is an Nx3
% matrix, x is an Mx3 matrix, inds is an Nx2 matrix specifying
% which azimuth and elevation elements a target point corresponds
% to.  azi and eli are vectors of length M.

chunk_size = 5000;

% Normalize x to be unit length
x = x ./ repmat(sqrt(sum(x.^2, 2)), 1, size(x,2));

% Rotate the sources by the opposite of the head rotation
R = euler_matrix(head_angles);
x = x * R';

% Find largest inner product with a target for each source.
N = size(x,1);
azi = zeros(N,1);
eli = zeros(N,1);

for i=1:chunk_size:N
  last = min(N, i+chunk_size-1);
  ip = target * x(i:last,:)';
  best = argmax(ip, 1);
  azi(i:last) = inds(best,1);
  eli(i:last) = inds(best,2);
end


function [h, cache] = lazy_hrtf(az, el, cache)
% Return the length-N FFT of the head-related impulse response at the
% specified azimuth and elevation.  Calculate as necessary.  Cache
% is a cell array indexed by {az,el} which contains impulse
% responses or their fourier transforms.

if isempty(cache.tfs{az,el})
  % Store resampling filter in cache so it doesn't need to be
  % recomputed every time, which is surprisingly expensive
  [resampled_ir, cache.firls_order] = ...
      resample(cache.irs{az,el}, cache.sr, 44100, cache.firls_order);
  % TODO: extra call to single()?
  cache.tfs{az,el} = single(fft(single(resampled_ir), cache.N));
end
h = cache.tfs{az,el};



function [head_width, els, azs, hrtf_cache, N, time_margin, pos3d, inds] = ...
    load_hrtf_data(cipicfile, N, sr)
% Load CIPIC hrir data from its file and initialize the aproporiate
% data structures.  hrtf_cache is a structure with fields irs,
% tfs, N, sr, firls_order.  Adds time_margin to N so that it can be
% trimmed off at the end of the calculation.

head_width = 0.145;
els        = [-45 : 5.625 : 230.625] * pi/180;
azs        = [-80 -65 -55 -45:5:45 55 65 80] * pi/180;

[A,E]   = meshgrid(azs, els);
[Ai,Ei] = meshgrid(1:length(azs), 1:length(els));
A = A(:); E = E(:); Ai = Ai(:); Ei = Ei(:);
pos3d = single([cos(A).*cos(E) sin(A) cos(A).*sin(E)]);
inds  = [Ai Ei];

% Load cipic hrtf, take length N fft.  Could be done outside and
% passed in
cipic = load(cipicfile);
hrir_l = permute(cipic.hrir_l, [3 1 2]);
hrir_r = permute(cipic.hrir_r, [3 1 2]);

% Pack hrirs into cells of hrtf cache
hrtf_cache.irs = cell(length(azs), length(els));
for azi = 1:length(azs)
  for eli = 1:length(els)
    hrtf_cache.irs{azi,eli} = [hrir_l(:,azi,eli) hrir_r(:,azi,eli)];
  end
end

% Remember how long the HRIRs are, and don't add any echos that are
% within that time of the end, they will wrap around to the
% beginning
time_margin = floor(length(hrtf_cache.irs{1,1}) * sr/44100);
N = N + time_margin;

% Initialize fields of hrtf_cache
hrtf_cache.tfs         = cell(length(azs), length(els));
hrtf_cache.N           = N;
hrtf_cache.sr          = sr;
hrtf_cache.firls_order = 10;



function [absorp, air_absorp] = ...
    setup_absorp(absorp, N, sr, temperature, humidity)
% Setup the absorption coefficients to either be 6x1 or 6xN.  If
% frequency-dependent, then interpolate at all frequencies

[r,c] = size(absorp);

% Interpolate frequencies for frequency-dependent absorption
if(c > 1)
  % Convert frequencies to [0,2pi]
  F = [125 250 500 1000 2000 4000]/sr * 2*pi;

  % Make the filters 0 at 0 and nyquist
  F = [0 F pi];
  z = ones(size(absorp,1),1);
  absorp = [z absorp absorp(:,end)];
  
  % interpolate the frequencies up to nyquist
  absorp = interp1(F, absorp', [0:N/2]'/N * 2*pi, 'cubic')';
  
  % reflect about nyquist to get the full spectrum
  absorp = absorp(:,[1:end end-1:-1:2]);

%   % Visualize absorption
%   subplot 211
%   plot([0:N-1]/N * sr, absorp')
%   subplot 212
%   plot(fftshift(ifft(absorp')))
%   drawnow
%   subplot 111
end

if r == 1
  absorp = repmat(absorp, 6, 1);
end

% Frequency-dependent absorption per meter traveled in the atmosphere.
% Gets multiplied by d, distance of source in meters
F_all = [0:N/2 N/2-1:-1:1]' * sr/N;
air_absorp = atmospheric_attenuation(1,temperature,humidity,F_all);
% Convert from dB to natural log
air_absorp = log(10)/20 * air_absorp;

% Cast to singles to speed up calculations a bit
absorp = single(absorp);
air_absorp = single(air_absorp);
