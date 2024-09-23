clc
clear
close all

% Simulation of an Arc-SAR system based on FMCW chirp

%% Anonymous Functions and constants


% carrier frequency
fc = 79e9;

c = physconst('LightSpeed');
lambda = c/fc;

% compute bw based on the range resolution
computeBw = @(range_res) c./(2.*range_res);

% compute the range based on the time of flight
range2time = @(range) (2.*range)./c; 

% compute the beat frequency based on the distance
range2beat = @(range, sweep) (sweep.*((2.*range))./c);

% compute azimuth resolution
az_res = @(r, beam_width) rad2deg(lambda./(2.*r.*deg2rad(beam_width)));


%% Init Radar Parameters

% I'm setting the ROI dimensions
ROI_l = 5;
ROI_w = 5;
% The scanned angle by the CSAR system
scan_angle = 180;

% I am assuming the maximum range is the same as the length of the ROI
range_max = ROI_l;
% the upsweept time depends on the maximum unambiguous range
% usually it is taken to be either 5 or 6 times greater the nominal value
% unlike the pulsed chirp, FMCW chirp does not have a PRI time
% it is not a monostatic system, we are using two antennas 
tm = 6*range2time(range_max); % this is correlated with the SNR and max dist
% tm = 68.9*10^-6;

% prf = 1e4;
% pri = 1./prf;
pri = tm;

range_res = 0.043;
ang_res = 0.5;
bw = computeBw(range_res);
sweep_slope = bw/tm;

% the sampling frequency of the A/D converter can be actually smaller than 
% twice the bandwidth. This holds true because of two reasons
% 1) complex sampled signal: the sample rate can be set the same as the bw
% 2) the sample rate can be computed based onto the max beat frequency 
% we should consider the sum of the range and doppler beat frequency,
% however the doppler beat frequency can be looked over since the targets
% are still
fb_max = range2beat(range_max, sweep_slope);
fs = max(2.*fb_max,bw);
fs = round(fs*tm)./tm;

%% FMCW chirp signal Model

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, ...
    'SampleRate',fs, 'SweepDirection', 'up');

figure
sig = waveform();
subplot(211); plot(0:1/fs:tm-1/fs,real(sig));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');

%% Radar Scene Model 

radius_platform = 0.13;

T = ((scan_angle/ang_res)*pri);
tStep = 0:pri:T;
t = 0:size(tStep, 2)-1;

% waypoints model for a circular trajectory
% initial positions of the radar platform
x0 = radius_platform;
y0 = 0;
z0 = 0;

% angular velocity for a semi-circle 
w = (pi)./T; 
% trajectory equation of a uniform circular motion
traj = w*tStep;  % theta = theta0 + w*t [radians]

x = radius_platform.*cos(traj);
y = radius_platform.*sin(traj);
z = z0*ones(size(t));
waypoints = [tStep.' x.' y.' z.'];

radarPlatform = phased.Platform('MotionModel','Custom','CustomTrajectory',waypoints);


%% Targets Model

num_targets = 3;
targetPos= ([1,2,0; -2,2,0; 2, 3, 0])';  
targetVel = zeros(size(targetPos));
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', ones(num_targets, 1)');
pointTargets = phased.Platform('InitialPosition', targetPos,'Velocity',targetVel);

%% Plot Ground-Truth

figure
plot(targetPos(1,1),targetPos(2,1),'ob')
hold on
plot(targetPos(1,2),targetPos(2,2),'ob')
hold on
plot(targetPos(1,3),targetPos(2,3),'ob')
hold on

% plot antenna beam pattern
th = linspace( 0, pi, 100);
R = range_max;
x_bp = R*cos(th);
y_bp = R*sin(th);
plot(x_bp,y_bp, '--r'); 
hold on

grid
axis equal
axis([-6 6 -0.5 6.5])
line1 = animatedline('DisplayName','Trajectory 1','Color','r','Marker','.');
title('Ground Truth')
addpoints(line1,x(1),y(1));

% Plot ROI
x_roi = [-5, 5, 5, -5, -5];
y_roi = [0.20, 0.20, 5, 5, 0.2];
plot(x_roi, y_roi, 'g-', 'LineWidth', 2);
hold on


% Plot turning table
for i = 1:length(t)
    step(radarPlatform, pri);
    addpoints(line1,x(i),y(i));
    pause(5*pri);
end

hold off;

radarPlatform.reset;


%% Radar Initialization 

antenna = phased.CosineAntennaElement('FrequencyRange', [79e9 83e9], 'CosinePower', [2, 2]);
az = -180:180;
el = -60:60;
pat = zeros(numel(el),numel(az),'like',1);
for m = 1:numel(el)
    temp = antenna(fc,[az;el(m)*ones(1,numel(az))]);
    pat(m,:) = temp;
end

%------------------- Visual-Test Rotating Antenna -------------------
figure;
pattern(antenna, fc, 'type', 'power')
index = 50;
newax = rotz(ang_res*(index-1));
rpat = rotpat(pat,az,el,newax);

txantenna = phased.CustomAntennaElement( ...
    'AzimuthAngles',az,'ElevationAngles',el,'SpecifyPolarizationPattern',false, ...
    'MagnitudePattern',mag2db(abs(rpat)), ...
    'PhasePattern',zeros(size(rpat)));

figure;
pattern(txantenna, fc, 'type', 'power')

% ----------------------------------------------------------------------

figure
imagesc(az,el,abs(pat))
axis xy
axis equal
axis tight
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
title('Original Radiation Pattern')
colorbar


%% ARC-SAR Simulation 

nSweep = length(t);
s = complex(zeros(waveform.SampleRate*waveform.SweepTime,nSweep));
s_unchirped = zeros(size(s));

for m = 1:nSweep

    [transmitter, radiator, channel, collector, receiver] = dynamicAntenna(fc, fs, pat, el, az, ang_res, m);

    % Update radar and target positions
    [radarPos,radarVel] = radarPlatform(pri);

%     restart(scene);
%     restart(radarPlat);
%     advance(scene);
%     pos = pose(radarPlat);
    
    [targetPos, targetVel] = pointTargets(pri);

    % Get the range and angle to the point targets
    [targetRange, targetAngle] = rangeangle(targetPos, radarPos);

    % Transmit FMCW waveform
    sig = waveform();
    txSig = transmitter(sig);

    % Radiate the pulse towards the targets
    txSig = radiator(txSig, targetAngle);
    
    % Propagate the signal and reflect off the target
    txSig = channel(txSig,radarPos,targetPos,radarVel,targetVel);

    % Reflect the pulse off the targets
    txSig = target(txSig);

    % Collect the reflected pulses at the antenna
    rxSig = collector(txSig, targetAngle);
    
    % Dechirp the received radar return
    rxSig = receiver(rxSig);    
    dechirpedSig = dechirp(rxSig,sig);
    
    % Visualize the spectrum
    s_unchirped(:,m) = rxSig;
    s(:,m) = dechirpedSig;
end

figure
imagesc(abs(s));
title('SAR Raw Data Dechirped')


%% Signal Sampling 

fs_adc = fb_max;
step_cut=round(size(s,1)/round(tm*fs_adc));
% s_sampled=s(1:step_cut:end,:);

%% Decimation 
% 
% Dn = fix(fs/(2*fb_max));
for m = size(s,2):-1:1
    s_dec(:,m) = decimate(s(:,m),step_cut,'FIR');
end

%% Backprojection Imaging

%fs = (2652e3); %ksps
ts = 1./fs_adc;
ts_start = 0;

i = floor(ROI_l./(range_res));
l = floor(scan_angle./(ang_res));

S = s_dec;

fastTime_sample = size(S,1);
slowTime_sample = size(S,2);
phi_m = linspace(0, scan_angle, slowTime_sample);
i_k = 1:fastTime_sample; 

A = @(R_t, phi_t) exp(1i*2*pi*(fc + sweep_slope*(i_k'-1)*ts+ts_start)*(2./c)*sqrt(R_t.^2 + radius_platform.^2 - 2*radius_platform*R_t*cosd(phi_m - phi_t)));

% image matrix ixl
g = zeros(i,l);

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

%% from Polar to Cartesian

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
%% Range Compression with matched filter

% correlation
[m, lg] = xcorr(s_unchirped(:,1), sig);
m = abs(m(lg > 0));
lm = lg(lg > 0);
[m1, lg1] = xcorr(s_unchirped(:,2), sig);
m1 = abs(m1(lg1 > 0));
lm1 = lg1(lg1 > 0);

% figure;
% subplot(2,2,1)
% plot(lm, m);
% title("First Pulse")
% subplot(2,2,2)
% plot(lm1, m1);
% title("Last Pulse")
% sgtitle("Range Compression with matched filter");

% convolution
s_fft = fft(s_unchirped, [], 1);
sig_fft = fft(conj(sig), [], 1);
sig_fft = repmat(sig_fft, [1 size(s_fft, 2)]);

rangeC_freq = sig_fft.*s_fft;

rangeC_time = abs(ifft(rangeC_freq, [], 1));

N = size(rangeC_time, 1);
faxis = linspace(0, tm, N);

% subplot(2,2,3)
% plot(faxis, rangeC_time(:,1));
% subplot(2,2,4)
% plot(faxis, rangeC_time(:,2));
% sgtitle("Range compression through convolution")

figure
imagesc((rangeC_time));
ylim([0 150]);


%% Beamwidth

figure;
% I'm printing the azimuth beamwidth with an elevation cut angle equal to 0
[beamWidth, angles] = beamwidth(antenna, fc, 'cut', 'Azimuth', 'Cutangle', 0);
beamwidth(antenna, fc, 'cut', 'Azimuth', 'Cutangle', 0)

SNR = radareqsnr(lambda, [5 5], 10, tm, 'Gain', [42 42]);
max_range = radareqrng(lambda,SNR,10,tm, 'unitstr', 'm');


%% --------------------- TESTS -------------------------

%% FFT
% figure;
% y = abs(fft(s_dec,[],1));
% N = size(y, 1);
% y = fftshift(y, 1);
% y = y/N;
% faxis = linspace(-fb_max/2, fb_max/2-1, N);
% faxis = beat2range(faxis', sweep_slope);
% figure
% plot(faxis, y(:,1), 'g')
% hold on
% plot(faxis, y(:,end), 'b')
% hold off;

