clc
clear
close all

% FMCW Simulation of a Circular SAR system

%% Anonymous Functions and constants

c = physconst('LightSpeed');

% compute bw based on the range resolution
computeBw = @(range_res) c./(2.*range_res);

% compute the range based on the time of flight
range2time = @(range) (2.*range)./c; 

% compute the beat frequency based on the distance
range2beat = @(range, sweep) sweep.*((2.*range)./c);


%% Init Parameters

% I'm setting the ROI dimensions
ROI_l = 5.5;
ROI_w = 180;

% carrier frequency
fc = 79e9;
lambda = c/fc;

% I am assuming the maximum range is the same as the width
range_max = ROI_l;
% the upsweept time depends on the maximum unambiguous range
% usually it is taken to be either 5 or 6 times greater the nominal value
% unlike the pulsed chirp, FMCW chirp does not have a PRI time
% it is not a monostatic system, we are using two antennas 
tm = 5.5*range2time(range_max);


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

%% FMCW chirp signal

waveform = phased.FMCWWaveform('SweepTime',tm,'SweepBandwidth',bw, ...
    'SampleRate',fs, 'SweepDirection', 'up');

sig = waveform();
subplot(211); plot(0:1/fs:tm-1/fs,real(sig));
xlabel('Time (s)'); ylabel('Amplitude (v)');
title('FMCW signal'); axis tight;
subplot(212); spectrogram(sig,32,16,32,fs,'yaxis');
title('FMCW signal spectrogram');

%% Radar Scene model 

format long

r_platform = 13e-2;

T = ((ROI_w/ang_res)*tm);
tStep = 0:tm:T;
t = 0:size(tStep, 2)-1;
vel = (2.*pi)./T;

% waypoints model for a circular trajectory
% initial positions of the radar platform
x0 = r_platform;
y0 = 0;
z0 = 0;

% angular velocity for a semi-circle 
angVel = (pi)./T; 
% trajectory equation of a uniform circular motion
law = angVel*tStep;  % theta = theta0 + w*t [radians]

x = r_platform.*cos(law);
y = r_platform.*sin(law);
z = z0*ones(size(t));
wpts = [tStep.' x.' y.' z.'];

radarPlatform = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts);

% figure;
% X = [];
% tstep = tm;
% nsteps = 361;
% for k = 1:nsteps
%     [pos,vel] = radarPlatform(tstep);
%     X = [X;pos'];    
% end
% plot(x,y,'o'); hold on
% plot(X(:,1),X(:,2),'.')
% hold off;

%% Target Model

num_targets = 3;
targetPos= ([1,2,0; -2,2,0; 4,3,0])'; 
targetVel = zeros(num_targets);
target = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', db2pow(10)*ones(num_targets, 1)');
pointTargets = phased.Platform('InitialPosition', targetPos,'Velocity',targetVel);

figure
plot(targetPos(1,1),targetPos(2,1),'*g')
hold on
plot(targetPos(1,2),targetPos(2,2),'*r')
hold on
plot(targetPos(1,3),targetPos(2,3),'*b')
hold on

% plot(X(:,1),X(:,2),'o');
% title('Configurazione SCENA')
xlim([-6 6])
ylim([-6 6])

%% Antenna parameters

% MODIFICARE F RANGE
txantenna =phased.IsotropicAntennaElement('FrequencyRange', [79e9 81e9]);
%txantenna = phased.CosineAntennaElement('FrequencyRange', [77e9 81e9]);
rxantenna = clone(txantenna);

% figure;
% pattern(txantenna,fc)

ant_aperture = 6.06e-4;                         % in square meter
ant_gain = aperture2gain(ant_aperture,lambda);  % in dB

tx_ppower = db2pow(5)*1e-3;                     % in watts
tx_gain = 9+ant_gain;                           % in dB

rx_gain = 15+ant_gain;                          % in dB
rx_nf = 4.5;                                    % in dB


peakpower = 10;
txgain = 36.0;
transmitter = phased.Transmitter( ...
    'PeakPower',tx_ppower, ...
    'Gain',tx_gain, ...
    'InUseOutputPort',true);
radiator = phased.Radiator( ...
    'Sensor',txantenna, ...
    'PropagationSpeed',c, ...
    'OperatingFrequency',fc);

channel = phased.FreeSpace( ...
    'SampleRate',fs, ...    
    'PropagationSpeed',c, ...
    'OperatingFrequency',fc, ...
    'TwoWayPropagation',true);

collector = phased.Collector( ...
    'Sensor',rxantenna, ...
    'PropagationSpeed',c, ...
    'OperatingFrequency',fc);
rxgain = 42.0;
noisefig = 10;
receiver = phased.ReceiverPreamp( ...
    'SampleRate',fs, ...
    'Gain',rx_gain, ...
    'NoiseFigure',rx_nf);

% specanalyzer = dsp.SpectrumAnalyzer('SampleRate',fs,...
%     'PlotAsTwoSidedSpectrum',true,...
%     'Title','Spectrum for received and dechirped signal',...
%     'ShowLegend',true);

%% CSAR Simulation 

nSweep = length(t);
s = complex(zeros(waveform.SampleRate*waveform.SweepTime,nSweep));
s_unchirped = zeros(size(s));

poss = [];

for m = 1:nSweep
    % Update radar and target positions
    [radarPos,radarVel] = radarPlatform(tm);
    
    poss = cat(1, poss, reshape(radarPos, [1, 3]));

    % THIS NEEDS TO BE FIXED
    % radarPos = X(m,:)';
    % radarVel = X_vel(m,:)';

    [targetPos, targetVel] = pointTargets(waveform.SweepTime);

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
    %specanalyzer([rxSig dechirpedSig]);
    s_unchirped(:,m) = rxSig;
    s(:,m) = dechirpedSig;
end

figure
imagesc(abs(s));
title('SAR Raw Data Dechirped')

%% Range Migration Correction for focusing

% vel_radarPlat = (pi.*r_platform)./(ROI_w/ang_res)*tm;
% vel = (2.*pi)./(ROI_w/ang_res)*tm;
% 
% slcimg = rangeMigrationFMCW(s,waveform,fc, vel_radarPlat, r_platform);
% 
% figure;
% imagesc(abs(slcimg));

%% Decimation 

% Dn = fix(fs/(2*fb_max));
% for m = size(s,2):-1:1
%     s_dec(:,m) = decimate(s(:,m),Dn,'FIR');
% end
% fs_d = fs/Dn;

%% Signal Sampling 

fs_adc= fb_max;
step_cut=round(size(s,1)/round(tm*fs_adc));
s_sampled=s(1:step_cut:end,:);


%% Range and Doppler estimation 

% rngdopresp = phased.RangeDopplerResponse('PropagationSpeed',c,...
%     'DopplerOutput','Speed','OperatingFrequency',fc,'SampleRate',fs,...
%     'RangeMethod','FFT','SweepSlope',sweep_slope,...
%     'RangeFFTLengthSource','Property','RangeFFTLength',2048,...
%     'DopplerFFTLengthSource','Property','DopplerFFTLength',256);
% 
% figure;
% plotResponse(rngdopresp,s_sampled);                     % Plot range Doppler map
% axis([-v_max v_max 0 range_max])

%% FFT tests

% s_fft = fft(s_dec, [], 1);
% s_fft = abs(s_fft);
% s_fft = s_fft./max(s_fft);
% % s_fft1 = beat2range(s_fft(:,1), sweep_slope);
% % s_fft2 = beat2range(s_fft(:,2), sweep_slope);
% 
% f = 1:size(s_fft, 1);
% range_raw=c*f/(2*sweep_slope);
% range=(ROI_l/range_raw(end))*range_raw';
% 
% figure;
% plot(range, s_fft(:,1), 'g'); 
% axis([0 5.5 0 1]);


%% FFT-RANGE
% mix_cut = slcimg;
% Mix_fft_2S=(fft(mix_cut,[],1)); %fft in matrice opera lungo colonne(X,[],1) per farlo lungo righe deve essere(X,[],2)
% Mix_fft_2S = fftshift((Mix_fft_2S/length(Mix_fft_2S)),1);%fft in matrice opera lungo colonne(X,1) per farlo lungo righe deve essere(X,2
% lm_2S=size(Mix_fft_2S,1);
% f_2S= fs_adc/lm_2S *(-lm_2S/2:lm_2S/2-1);
% figure
% plot(abs(Mix_fft_2S(:,1)));
% hold on
% plot(abs(Mix_fft_2S(:,end)),'r');
% title('FFT-DS')
% % 
% % FFT Singole Side
% lm=round(lm_2S/2);
% Mix_fft = 2*Mix_fft_2S(lm-1:-1:1,:);
% f=(f_2S(round(lm_2S/2)+1:end)); 
% range_raw=c*f/(2*sweep_slope);
% range=(maxRange/range_raw(end))*range_raw';
% 
% figure
% plot(range(1:end),abs(Mix_fft(:,1)));
% hold on
% plot(range(1:end),abs(Mix_fft(:,end)),'r');
% % plot(abs(Mix_fft(1:end,:)))
% title('FFT-SS')
% xlim([0 6])
% ylabel('Amplitude')

%% Imaging through Back-Projection 1

%    l
%   ---------
% i |            Number of Pixels of the image
%   |              
i = floor(ROI_l./(range_res));
l = floor(ROI_w./(ang_res)+1);


S = s_sampled;
fastTime_sample = size(S,1);

phi_m = 0:ang_res:180;
i_k = 0:fastTime_sample-1; % (i-1)


% tau = 2*R/c -> R is the distance from the radar to the target 
% R is a row 1xM, where M is the number of slow-time samples
tau = @(r_l, phi_l, phi_m) sqrt(r_l.^2+r_platform.^2-2.*r_platform.*r_l*cosd(phi_m-phi_l)).*(2./c);
a = @(r_l, phi_l, phi_m, i_k) exp(1i.*2.*pi.*(fc.*tau(r_l, phi_l, phi_m) - (sweep_slope./2)*tau(r_l, phi_l, phi_m).^2 + i_k'*(bw./fastTime_sample)*tau(r_l, phi_l, phi_m)));

% image matrix ixl
g = zeros(i,l);

tic 
for q = 1:l
    for p = 1:i
        a_il_cs = a((p-1)*range_res, (q-1)*ang_res, phi_m, i_k);
        g(p,q) = conj(reshape(a_il_cs,[1,numel(a_il_cs)]))*reshape(S, [numel(S), 1]);
     end
end
toc

figure;
imagesc(abs(g));

%% from Polar to Cartesian

r_l = (0:range_res:(ROI_l-range_res))';  % Create column vector for range
theta_l = 0:ang_res:ROI_w;  % Create row vector for angles
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
xq = linspace(min(Xc_flat), max(Xc_flat), 128);
yq = linspace(min(Yc_flat), max(Yc_flat), 361);
[Xq, Yq] = meshgrid(xq, yq);
% Interpolate the intensity values on the Cartesian grid
Interpolated_intensity = Cartesian_intensity(Xq, Yq);
% Plot the result (optional)
figure;
imagesc(xq, yq, abs(Interpolated_intensity));
axis xy;
xlabel('X (meters)');
ylabel('Y (meters)');
title('Interpolated Intensity in Cartesian Coordinates');
colorbar;