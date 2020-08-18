clear;
%% Get Files
% Get fMax
fileID_fMax = fopen('fMax.txt', 'r');
formatSpec = '%f';
fMax = fscanf(fileID_fMax, formatSpec);
% Open files
loadFile = strcat("load", num2str(fMax), ".txt");
nodeFile = strcat("nodeOut", num2str(fMax), ".txt");
timeFile = strcat("time", num2str(fMax), ".txt");
sigmaFile = strcat("sigma", num2str(fMax), ".txt");
fileID_load = fopen(loadFile,'r');
fileID_nodeOut = fopen(nodeFile,'r');
fileID_time = fopen(timeFile,'r');
fileID_sigma = fopen(sigmaFile,'r');
load = fscanf(fileID_load, formatSpec);
node = fscanf(fileID_nodeOut, formatSpec);
t = fscanf(fileID_time, formatSpec);
t0 = t(1:length(load)-1);
sigma = fscanf(fileID_sigma, formatSpec);
fclose(fileID_load);
fclose(fileID_nodeOut);
fclose(fileID_time);
fclose(fileID_sigma);
%% Even Length Signals
if rem(length(load), 2) == 1
load(end) = [];
node(end) = [];
t(end) = [];
end
%% Filter Signal
passband = 20;
filtered = highpass(node, passband, length(node));
max_gain = 0.1;
scaled = rescale(filtered, max_gain/max(filtered)*min(filtered), max_gain);
%% Response FFT
Y = fft(node);
Fs = length(node)/2;
T = 1/Fs;
L = floor(length(node)/2)*2;
P2 = abs(Y/L);
f0 = Fs*(0:(L-1))/L;
P1 = P2(1:L/2);
P1(2:end-1) = 2*P1(2:end-1);
f = linspace(0, Fs, length(P1));
%% Convolute IR and New Signal
[data,fs_new]=audioread('Master_Of_Puppets.wav');
w = conv(data(1:500000,1), scaled, 'valid');
%% Plot
figure(1); clf
% Load
subplot(2, 1, 1)
plot(t0, load)
title('Load')
grid on
grid minor
xlabel('Time (s)')
ylabel('Pressure (Pa)')
% Output
subplot(2, 1, 2)
plot(t, node)
title('Mic Response')
grid on
grid minor
xlabel('Time (s)')
ylabel('Pressure (Pa)')
% FFT
figure(2); clf
subplot(2, 1, 1)
plot(f,P1)
grid on
title('Single-Sided Amplitude Spectrum of Node(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
subplot(2, 1, 2)
plot(t, scaled)
grid on
title('Impulse Response with Frequencies lower than 20Hz Omitted')
xlabel('Time (s)')
ylabel('Pressure (Pa)')
