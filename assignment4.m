% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 04
% Analysis, Modeling, and Resynthesis of Music Signals
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script performs LPC-based analysis and resynthesis of a
% guitar signal using two methods:
%   1) Noise-excited LPC
%   2) Residual-excited LPC
% ===============================================================

clc;
clear;
close all;

%% -------------------- Load Audio --------------------------
[fileData, Fs] = audioread('guitar_track.mp3');

% Convert stereo to mono if needed
if size(fileData, 2) > 1
    x = mean(fileData, 2);
else
    x = fileData;
end

% Normalize signal
x = x / max(abs(x));

%% -------------------- Select Segment ----------------------
% Use a short clean segment for analysis
start_sec = 0;          % change if needed
dur_sec   = 8;          % 5 to 10 seconds is enough

start_sample = floor(start_sec * Fs) + 1;
end_sample   = min(length(x), floor((start_sec + dur_sec) * Fs));

x = x(start_sample:end_sample);

% Normalize again after trimming
x = x / max(abs(x));

t = (0:length(x)-1) / Fs;

%% -------------------- LPC Analysis ------------------------
% LPC order can be adjusted depending on the signal
p = 24;

% Estimate LPC coefficients
a = lpc(x, p);

%% -------------------- Method 1: Noise-Excited LPC ---------
% White noise is used as the excitation signal
excitation_noise = randn(length(x), 1);
excitation_noise = excitation_noise / max(abs(excitation_noise));

y_noise = filter(1, a, excitation_noise);
y_noise = y_noise / max(abs(y_noise));

%% -------------------- Method 2: Residual-Excited LPC ------
% Extract the residual from the original signal
residual = filter(a, 1, x);

% Re-synthesize using the residual
y_residual = filter(1, a, residual);
y_residual = y_residual / max(abs(y_residual));

%% -------------------- Plot 1: Time-Domain Comparison ------
fig1 = figure;

subplot(3,1,1);
plot(t, x, 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Guitar Signal');
grid on;

subplot(3,1,2);
plot(t, y_noise, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Noise-Excited LPC Output');
grid on;

subplot(3,1,3);
plot(t, y_residual, 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('Residual-Excited LPC Output');
grid on;

saveas(fig1, 'guitar_lpc_time_comparison.png');

%% -------------------- Plot 2: Overlay Comparison ----------
fig2 = figure;

plot(t, x, 'b', 'LineWidth', 1); hold on;
plot(t, y_noise, 'r');
plot(t, y_residual, 'g');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs LPC Resynthesized Outputs');
legend('Original', 'Noise LPC', 'Residual LPC');
grid on;

saveas(fig2, 'guitar_lpc_overlay.png');

%% -------------------- Plot 3: Frequency-Domain Comparison -
N = length(x);
f = (0:N-1) * (Fs / N);

X_original = abs(fft(x));
X_noise    = abs(fft(y_noise));
X_residual = abs(fft(y_residual));

fig3 = figure;

plot(f, X_original, 'b'); hold on;
plot(f, X_noise, 'r');
plot(f, X_residual, 'g');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Frequency-Domain Comparison');
legend('Original', 'Noise LPC', 'Residual LPC');
grid on;

saveas(fig3, 'guitar_lpc_frequency_comparison.png');

%% -------------------- Save Audio Outputs ------------------
audiowrite('guitar_original_segment.wav', x, Fs);
audiowrite('guitar_noise_lpc.wav', y_noise, Fs);
audiowrite('guitar_residual_lpc.wav', y_residual, Fs);

%% -------------------- Play Results ------------------------
disp('Playing original segment...');
soundsc(x, Fs);
pause(length(x)/Fs + 1);

disp('Playing noise-excited LPC output...');
soundsc(y_noise, Fs);
pause(length(y_noise)/Fs + 1);

disp('Playing residual-excited LPC output...');
soundsc(y_residual, Fs);
