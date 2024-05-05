clear;
clc;
close all;

%% Anonymous functions and Constants 

c = physconst('LightSpeed');

% PRI is computed according to the maximum unambiguous range equation
comp_PRI = @(d_max) (2.*d_max)./c;

% Maximum distance that can be scanned by the radar - it depends on fs
comp_dMax = @(fs, b, tp) (fs.*c)./(2.*(b./tp));

% PulseWidth of the chirp, it is the upsweep time of the chirp
comp_Tp = @(dmax, B, fs) (2.*dmax.*B)./(fs.*c);

% Range resolution - it depends on the bw
comp_rangeRes = @(b) c./(2.*b);

%% Initialize Radar and Environment Parameters

% I'm setting the ROI dimensions - ROI_l is the maximum distance
ROI_w = 900;
ROI_l = 900;
% radius of the physical platform the radar is mounted on
r_platform = 0.13;

% fast time sampling frequency
fs = 150e6;
% center(carrier) frequency
fc = 77e9;
pri = 7e-6;
% slow time sampling frequency
prf = 1/pri;
% pulsewidth of the upsweeft chirp
tp =  1.4e-7;

bw = fs./2;

% Range and Angle Resolution
range_res = comp_rangeRes(bw);
ang_res = 0.45;

%% LFM Pulse Generator
% A LFM pulse wave is generated
% #tot_samples -> Fs/PRF = PRI/Ts -> we don't need all of these
% #pulse_samples -> PRI/TS - (PRI-PW)/Ts -> these are the ones we need

dutyCycle = tp./pri;
wave = phased.LinearFMWaveform('SampleRate',fs, 'PulseWidth', tp, ...
    'PRF', prf, 'SweepBandwidth', bw, 'NumPulses',1);

sig = wave();
Nr = length(sig);

figure('Name','Chirp wave','NumberTitle','off');
wav = step(wave);
numpulses = size(wav,1);
t = (0:(numpulses-1))/wave.SampleRate;
plot(t*1e6,real(wav))
xlabel('Time (\mu sec)')
ylabel('Amplitude')

%% Scenario - Targets
% A radar scenario simulates a 3-D environment containing multiple platforms
% I'm considering a scenario composed of three different targets
% Such targets are located 500, 530, 750 meters away from the radar on the
% x - axis, their speed is supposed to be zero, hence they are fixed in
% space
% We consider targets whose cross section is 10dB


% SAR_sim = radarScenario;
% target_1 = platform(SAR_sim);
% target_2 = platform(SAR_sim);

Numtgts = 3;
tgtpos = zeros(Numtgts);
tgtpos(1,:) = [500 800 850];
tgtvel = zeros(3,Numtgts);
tgtvel(1,:) = [0 0 0];
tgtrcs = db2pow(10)*[1 1 1];
tgtmotion = phased.Platform(tgtpos,tgtvel);
target = phased.RadarTarget('PropagationSpeed',c,'OperatingFrequency',fc, ...
    'MeanRCS',tgtrcs);

figure('Name','Target Positions in the scene','NumberTitle','off');
h = axes;plot(tgtpos(2,1),tgtpos(1,1),'*g');hold all;plot(tgtpos(2,2),tgtpos(1,2),'*r');hold all;plot(tgtpos(2,3),tgtpos(1,3),'*b');hold off;
set(h,'Ydir','reverse');xlim([-ROI_l ROI_l]);ylim([0 ROI_l]);
title('Ground Truth');ylabel('Range');xlabel('Cross-Range');

%% Scenario - SAR

radarpos = [0;0;0];
radarvel = [0;0;0];
radarmotion = phased.Platform(radarpos,radarvel);
% Waypoints model
x0 = 0;
y0 = ROI_w;
z0 = 0;
vx = 5;
vy = 10;
vz = 0;
ax = 1;
ay = -1;

t = [0:1:400];
x = x0 + vx*t + ax/2*t.^2;
y = y0 + vy*t + ay/2*t.^2;
z = z0*ones(size(t));
wpts = [t.' x.' y.' z.'];

pltfm = phased.Platform('MotionModel','Custom','CustomTrajectory',wpts);
tstep = .5;
nsteps = 41;
X = [];

for k = 1:nsteps
    [pos,vel] = pltfm(tstep);
    X = [X;pos'];   
    figure;
    plot(x,y,'o'); hold on
    plot(X(:,1),X(:,2),'.')
    hold off;
end

%% Rx and Tx Antennas Front-end and Back-end
% They are both placed at the origin (monostatic)

txantenna = phased.IsotropicAntennaElement;
rxantenna = clone(txantenna);

peakpower = 10;
txgain = 36.0;
transmitter = phased.Transmitter( ...
    'PeakPower',peakpower, ...
    'Gain',txgain, ...
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
    'Gain',rxgain, ...
    'NoiseFigure',noisefig);

%% Data cube generation
% I'm generating a radar cube NrxNp
%
% Nr -> #samples of the chirp -> fast-time
% Np -> #pulses(step angles Tx is transmitting towards to) -> slow-time
% Pulses are being transmitted towards only one direction in this set-up
%
% We first compute the actual value of the PRI of a pulse that is needed in
% order to cover a predefined maximum unambiguous distance. Then the Nr, 
% which is the actual number of samples, is computed

%Nr = ceil(comp_PRI(ROI_l).*fs);
Np = 128;
radarCube = zeros(Nr,Np);
for n = 1:Np
    [sensorpos,sensorvel] = radarmotion(pri);
    [tgtpos,tgtvel] = tgtmotion(pri);
    [tgtrng,tgtang] = rangeangle(tgtpos,sensorpos);
    sig = wave();
    % we use only the pulse length that covers the targets
    sig = sig(1:Nr);
    [txsig,txstatus] = transmitter(sig);
    txsig = radiator(txsig,tgtang);
    txsig = channel(txsig,sensorpos,tgtpos,sensorvel,tgtvel);    
    tgtsig = target(txsig);   
    rxcol = collector(tgtsig,tgtang);
    rxsig = receiver(rxcol);
    radarCube(:,n) = rxsig;
end

%% Display the image of the data cube containing signals per pulse

figure;
imagesc((0:(Np-1))*pri*1e6,(0:(Nr-1))/fs*1e6,abs(radarCube))
xlabel('Slow Time {\mu}s')
ylabel('Fast Time {\mu}s')

%% Range Response - Range FFT
matchingcoeff = getMatchedFilter(wave);
ndop = 1;
rangeresp = phased.RangeResponse('SampleRate',fs,'PropagationSpeed',c);
[resp,rnggrid] = rangeresp(radarCube,matchingcoeff);
figure;
imagesc((1:Np),rnggrid,abs(resp))
xlabel('Pulse')
ylabel('Range (m)')

%% Integrate noncoherent
figure;
intpulse = pulsint(resp(:,1:20),'noncoherent');
plot(rnggrid,abs(intpulse))
xlabel('Range (m)')
title('Noncoherent Integration of 20 Pulses')

%% BackProjection

% bpa_processed = helperBackProjection(resp,rnggrid,ceil(Nr./fs),fc,fs,prf,0,(c./(2.*bw)),c);
% 
% figure()
% imagesc((abs(bpa_processed(600:1400,1700:2300))));
% title('SAR Data focused using Back-Projection algorithm ')
% xlabel('Cross-Range Samples')
% ylabel('Range Samples')