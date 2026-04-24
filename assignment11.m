% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 11
% Spectral Estimation
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Description:
% This script demonstrates the difference between ordinary FFT
% analysis and Welch spectral estimation using a synthetic chord.
%
% NOTE:
% - MATLAB Online friendly
% - figures are shown on screen
% - please take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  Signal generation
% --------------------------------------------------------------
fs = 4000;                 % Sampling frequency
t  = 0:1/fs:0.5;           % Short duration to worsen FFT resolution
t  = t(:);

% Frequencies intentionally chosen between FFT bins
f1 = 445.3;
f2 = 612.7;
f3 = 789.2;

x = sin(2*pi*f1*t) + sin(2*pi*f2*t) + sin(2*pi*f3*t);
x = x / max(abs(x) + eps);

N = length(x);

%% -------------------------------------------------------------
%  Reduced vectors for time-domain plotting
% --------------------------------------------------------------
step_t = max(1, floor(N / 2000));
idx_t  = 1:step_t:N;

t_plot = t(idx_t);
x_plot = x(idx_t);

%% -------------------------------------------------------------
%  FFT
% --------------------------------------------------------------
X = fft(x);

halfN = floor(N/2);
f_axis = (0:halfN-1) * (fs/N);
X_mag = abs(X(1:halfN));

%% -------------------------------------------------------------
%  Welch spectral estimation
% --------------------------------------------------------------
window = 256;
[pxx, w] = pwelch(x, window, [], [], fs);

%% -------------------------------------------------------------
%  Figure 1: Original chord signal in time domain
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 420]);
plot(t_plot, x_plot, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Chord Signal (Time Domain)');
legend('Chord Signal');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 2: FFT vs Welch spectral estimation
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 650]);

subplot(2,1,1);
plot(f_axis, X_mag, 'r', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum using FFT');
legend('FFT');
xlim([0 1200]);
grid on;

subplot(2,1,2);
plot(w, pxx, 'b', 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Spectrum using Welch Spectral Estimation');
legend('Welch');
xlim([0 1200]);
grid on;

drawnow;

%% -------------------------------------------------------------
%  Save audio output
% --------------------------------------------------------------
audiowrite('assignment11_synthetic_chord.wav', x, fs);

disp('Assignment 11 completed successfully.');
disp('Please take screenshots of the 2 figures manually.');
disp('Suggested screenshot names:');
disp('  assignment11_time_domain.png');
disp('  assignment11_fft_vs_welch.png');
disp('Audio file saved:');
disp('  assignment11_synthetic_chord.wav');
