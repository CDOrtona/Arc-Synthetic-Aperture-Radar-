clear
clc;
close all;


% 
% This algorithm is based onto the simulation algorithm but it uses real data acquired with the radar sensor AWR1642
% The imaging algorithm is based onto the Backprojection method
% There will be implemented two versions of Backprojection:
%     The first version of the algorithm has a problem with the side lobes, hence the resulting image has bda artifacts 
%     The second version implements a slight modification of the algorithm which helps to reduce almost completly all the artifacts 


%% Anonymous Functions and constants


% carrier frequency
fc = 77e9;

c = physconst('LightSpeed');
lambda = c/fc;

% compute bw based on the range resolution
computeBw = @(range_res) c./(2.*range_res);

% compute range res based on the bw
compute_rangeRes = @(bw) c./(2.*bw);

% compute the time of flight from range
range2time = @(range) (2.*range)./c; 

% compute range from time
time2range = @(time) (c*time)./2;

% compute the beat frequency based on the distance
range2beat = @(range, sweep) sweep.*((2.*range)./c);

% compute azimuth resolution
az_res = @(r, beam_width) rad2deg(lambda./(2.*r.*deg2rad(beam_width)));

% compute azimuth bandwidth
az_band = (4*0.45*sind(120/2))/lambda;

%% Raw Radar Data Read Function

S = dca1000ReadRaw("sim/Prova_Esterna_10.bin");


%% Init Parameters

% I'm setting the ROI dimensions
ROI_l = 10;
ROI_w = 10;
% The scanned angle by the Arc-SAR system
scan_angle = 180;

% I am assuming the maximum range is the same as the length of the ROI
range_max = ROI_l;
% the upsweept time depends on the maximum unambiguous range
% usually it is taken to be either 5 or 6 times greater the nominal value
% unlike the pulsed chirp, FMCW chirp does not have a PRI time
% it is not a monostatic system, we are using two antennas 
%tm = 6*range2time(range_max); % this is correlated with the SNR and max dist
% tm = 68.9*10^-6;

% Beamwidth of the antennas mounted on the Radar
hpbw = 60;
half_hpbw = hpbw./2;

% Radius of the physical 3D printed rotating platform
radius_platform = 0.20; 

%bandwidth
bw = 3999.54e6;

range_res = compute_rangeRes(bw);
ang_res = az_res(radius_platform, 65);

% Chirp rate
sweep_slope = 39.058e12;

% Time duration of the Chirp
tm  = bw/sweep_slope;
pri = tm;

% the sampling frequency of the A/D converter can be actually smaller than 
% twice the bandwidth. This holds true because of two reasons
% 1) complex sampled signal: the sample rate can be set the same as the bw
% 2) the sample rate can be computed based onto the max beat frequency 
% we should consider the sum of the range and doppler beat frequency,
% however the doppler beat frequency can be looked over since the targets
% are still
fb_max = range2beat(range_max, sweep_slope);
% fs = max(2.*fb_max,bw);
% fs = round(fs*tm)./tm;


%% Back Projection Imaging - This produces artifacts in the radar image

% sampling frequency 
fs = (2652e3); %ksps
ts = 1./fs;
ts_start = 3.94e-6;

% dimensions in pixels of the output radar image
i = floor(ROI_l./(range_res));
l = floor(scan_angle./(ang_res));

fastTime_sample = size(S,1);
slowTime_sample = size(S,2);
phi_m = linspace(0, scan_angle, slowTime_sample);
i_k = 1:fastTime_sample; 

% Anonymous Function: It computes the matrix A which is the NxM matrix
% where N is the fast-time samples and M is the slow-time samples
A = @(R_t, phi_t) exp(1i*2*pi*(fc + sweep_slope*(i_k'-1)*ts+ts_start)*(2./c)*sqrt(R_t.^2 + radius_platform.^2 - 2*radius_platform*R_t*cosd(phi_m - phi_t)));

% image matrix ixl
g = zeros(i,l);

% This loop computes all the ixl pixel of an image. In order to compute
% each one of them, a matrix A must be generated and the matched filter
% operation based on the Observation Model S = Ag must be computed
tic 
for q = 1:l
    for pos = 1:i
        if((pos-1)*range_res >= radius_platform)
            a_il = A(((pos-1)*range_res), ((q-1)*ang_res));
            g(pos,q) = conj(reshape(a_il,[1,numel(a_il)]))*reshape(S, [numel(S), 1]);
        else
            g(pos,q) = 0;
        end
     end
end
toc

figure;
imagesc(abs(g));
colormap jet;
colorbar;

%% Backprojection Imaging Narrower Aperture 

fs = (2652e3); %ksps
ts = 1./fs;
ts_start = 3.94e-6;

i = floor(ROI_l./(range_res));
l = floor(scan_angle./(ang_res));

fastTime_sample = size(S,1);
slowTime_sample = size(S,2);
phi_m = linspace(0, scan_angle, slowTime_sample);
i_k = 1:fastTime_sample; 

A = @(R_t, phi_t, phi_m) exp(1i*2*pi*(fc + sweep_slope*(i_k'-1)*ts+ts_start)*(2./c)*sqrt(R_t.^2 + radius_platform.^2 - 2*radius_platform*R_t*cosd(phi_m - phi_t)));

% image matrix ixl
g = zeros(i,l);

% this for loop does the same thing as the previous one with the only
% difference that for each pixel there is no need to use the whole matrix A
% but rather a more efficient and less time consuming matrix A' is used
% This new matrix A has a fewer number of fast-time samples which are
% computed dinamically according to the specific pixel which is being
% computed by the algorithm
tic 
for q = 1:l
    for pos = 1:i
        if((pos-1)*range_res >= 1)
            phi_v = phi_m - (q-1)*ang_res;
            n_ind = find(abs(phi_v) < half_hpbw);
            s_temp = S(:, n_ind);
            phi_m_temp = phi_m(n_ind);
            a_il = A(((pos-1)*range_res), ((q-1)*ang_res), phi_m_temp);
            g(pos,q) = conj(reshape(a_il,[1,numel(a_il)]))*reshape(s_temp, [numel(s_temp), 1]);
        else
            g(pos,q) = 0;
        end
     end
end
toc


figure;
imagesc(abs(g));
colormap jet;
colorbar;

%% from Polar to Cartesian coordinates conversion

r_l = (0:range_res:(ROI_l-range_res))';  % Create column vector for range
theta_l = 0:ang_res:scan_angle;  % Create row vector for angles
theta_l = theta_l(1:end-1);
% Convert polar to Cartesian coordinates
[Theta, R] = meshgrid(theta_l, r_l);
[Xc, Yc] = pol2cart(deg2rad(Theta), R);  % Use radians for trigonometric functions
% Ensure Intensity_polarMatrix is of the correct size
%assert(all(size(g) == size(Theta)), 'Size mismatch between Intensity_polarMatrix and Theta');
% Flatten the matrices for interpolation
Xc_flat = Xc(:);
Yc_flat = Yc(:);
Intensity_flat = g(:);
% Create scattered interpolant
Cartesian_intensity = scatteredInterpolant(Xc_flat, Yc_flat, Intensity_flat, 'linear', 'none');
% Create a Cartesian grid for interpolation
xq = linspace(min(Xc_flat), max(Xc_flat), length(r_l));
yq = linspace(min(Yc_flat), max(Yc_flat), length(theta_l));
[Xq, Yq] = meshgrid(xq, yq);
% Interpolate the intensity values on the Cartesian grid
Interpolated_intensity = Cartesian_intensity(Xq, Yq);
% Plot the result (optional)
figure;
imagesc(xq, yq, abs(Interpolated_intensity));
axis xy;
xlabel('X (meters)');
ylabel('Y (meters)');
title('Back-Projection Imaging');
colormap jet;
colorbar;
hold on;

th = linspace( 0, pi, 100);
for ii = 1:ROI_w
    R = 1;
    x = ii.*R*cos(th);
    y = ii*R*sin(th);
    plot(x, y, '--', 'Color', '#bf3918'); 
end

