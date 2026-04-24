% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 09
% Music Source Separation
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Description:
% This script demonstrates separation of two mixed musical
% sources using:
%   1) PCA-based method
%   2) ICA / matrix inversion method
%   3) Band-pass filter method
%
% NOTE:
% - Figures are shown on screen
% - Please take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  Basic parameters
% --------------------------------------------------------------
fs = 8000;                 % Sampling frequency
dur = 2;                   % Duration in seconds
t = 0:1/fs:dur;
t = t(:);

%% -------------------------------------------------------------
%  Create two source signals (two different chords)
% --------------------------------------------------------------
% Source 1: lower-frequency chord
f1 = [261.63, 329.63, 392.00];   % C major chord
S1 = sin(2*pi*f1(1)*t) + sin(2*pi*f1(2)*t) + sin(2*pi*f1(3)*t);

% Source 2: higher-frequency chord
f2 = [659.25, 783.99, 987.77];   % Higher octave chord
S2 = sin(2*pi*f2(1)*t) + sin(2*pi*f2(2)*t) + sin(2*pi*f2(3)*t);

% Normalize sources
S1 = S1 / max(abs(S1) + eps);
S2 = S2 / max(abs(S2) + eps);

%% -------------------------------------------------------------
%  Mix signals into two microphone observations
% --------------------------------------------------------------
A = [0.6 0.3;
     0.4 0.7];

Y = A * [S1.'; S2.'];   % 2 x N matrix
y1 = Y(1,:).';
y2 = Y(2,:).';

mix = y1 + y2;
mix = mix / max(abs(mix) + eps);

%% -------------------------------------------------------------
%  Prepare reduced time vectors for plotting
% --------------------------------------------------------------
N = length(mix);
step_t = max(1, floor(N / 2000));
idx_t  = 1:step_t:N;

t_plot   = t(idx_t);
y1_plot  = y1(idx_t);
y2_plot  = y2(idx_t);
mix_plot = mix(idx_t);
S1_plot  = S1(idx_t);
S2_plot  = S2(idx_t);

%% =============================================================
%  1) PCA-based separation
% =============================================================
X_pca = [y1 y2];
X_pca = X_pca - mean(X_pca, 1);

C = cov(X_pca);
[V,D] = eig(C);

[eigvals, idx] = sort(diag(D), 'descend');
V = V(:, idx);

Y_pca = X_pca * V;
PC1 = Y_pca(:,1);
PC2 = Y_pca(:,2);

PC1 = PC1 / max(abs(PC1) + eps);
PC2 = PC2 / max(abs(PC2) + eps);

variance_explained = eigvals / sum(eigvals + eps);

%% -------------------------------------------------------------
%  Figure 1: Mixed microphone signals
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 500]);

subplot(2,1,1);
plot(t_plot, y1_plot, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Microphone Signal y_1');
grid on;

subplot(2,1,2);
plot(t_plot, y2_plot, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Microphone Signal y_2');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 2: PCA variance explained
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 800 450]);
bar(variance_explained);
xlabel('Principal Component');
ylabel('Variance Ratio');
title('Variance Explained by Principal Components');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 3: PCA separated signals
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 450]);
plot(t_plot, PC1(idx_t), 'r', 'LineWidth', 1); hold on;
plot(t_plot, PC2(idx_t), 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Separated Signals using PCA');
legend('Component 1','Component 2');
grid on;

drawnow;

%% =============================================================
%  2) ICA / Matrix inversion based separation
% =============================================================
A_inv = inv(A);
S_est = A_inv * Y;

S_est1 = S_est(1,:).';
S_est2 = S_est(2,:).';

S_est1 = S_est1 / max(abs(S_est1) + eps);
S_est2 = S_est2 / max(abs(S_est2) + eps);

%% -------------------------------------------------------------
%  Figure 4: Mixing matrix
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 800 450]);
bar(A);
xlabel('Source Index');
ylabel('Mixing Coefficient');
title('Mixing Matrix A');
legend('Microphone 1', 'Microphone 2');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 5: ICA recovered signals
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 450]);
plot(t_plot, S_est1(idx_t), 'r', 'LineWidth', 1); hold on;
plot(t_plot, S_est2(idx_t), 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Signals using ICA / Matrix Inversion');
legend('Recovered Source 1','Recovered Source 2');
grid on;

drawnow;

%% =============================================================
%  3) Band-pass-filter-based separation
% =============================================================
[b1,a1] = butter(6, [200 450]/(fs/2), 'bandpass');
[b2,a2] = butter(6, [600 1100]/(fs/2), 'bandpass');

sep1 = filter(b1, a1, mix);
sep2 = filter(b2, a2, mix);

sep1 = sep1 / max(abs(sep1) + eps);
sep2 = sep2 / max(abs(sep2) + eps);

% FFT prep
MixFFT  = abs(fft(mix));
SepFFT1 = abs(fft(sep1));
SepFFT2 = abs(fft(sep2));

halfN = floor(N/2);
f_axis = (0:halfN-1) * (fs/N);

MixFFT_half  = MixFFT(1:halfN);
SepFFT1_half = SepFFT1(1:halfN);
SepFFT2_half = SepFFT2(1:halfN);

step_f = max(1, floor(halfN / 2500));
idx_f  = 1:step_f:halfN;

f_plot = f_axis(idx_f);
Mix_plot  = MixFFT_half(idx_f);
Sep1_plot = SepFFT1_half(idx_f);
Sep2_plot = SepFFT2_half(idx_f);

%% -------------------------------------------------------------
%  Figure 6: Mixed signal in time and frequency
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 500]);

subplot(2,1,1);
plot(t_plot, mix_plot, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Mixed Signal - Time Domain');
grid on;

subplot(2,1,2);
plot(f_plot, Mix_plot, 'LineWidth', 1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Mixed Signal - Frequency Domain');
xlim([0 1200]);
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 7: Filter responses and spectra
% --------------------------------------------------------------
[H1,w1] = freqz(b1, a1, 512, fs);
[H2,w2] = freqz(b2, a2, 512, fs);

figure('Color','w','Position',[100 100 1000 500]);

subplot(2,1,1);
plot(w1, abs(H1), 'r', 'LineWidth', 1.2); hold on;
plot(w2, abs(H2), 'b', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Band-Pass Filter Responses');
legend('Filter for Source 1', 'Filter for Source 2');
grid on;

subplot(2,1,2);
plot(f_plot, Sep1_plot, 'r', 'LineWidth', 1.0); hold on;
plot(f_plot, Sep2_plot, 'b', 'LineWidth', 1.0);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectra of Filter-Separated Signals');
legend('Separated Source 1','Separated Source 2');
xlim([0 1200]);
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 8: Band-pass separated signals
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 450]);
plot(t_plot, sep1(idx_t), 'r', 'LineWidth', 1); hold on;
plot(t_plot, sep2(idx_t), 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Separated Signals using Band-Pass Filters');
legend('Separated Source 1','Separated Source 2');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment9_source1.wav', S1, fs);
audiowrite('assignment9_source2.wav', S2, fs);
audiowrite('assignment9_mixed_signal.wav', mix, fs);

audiowrite('assignment9_pca_component1.wav', PC1, fs);
audiowrite('assignment9_pca_component2.wav', PC2, fs);

audiowrite('assignment9_ica_source1.wav', S_est1, fs);
audiowrite('assignment9_ica_source2.wav', S_est2, fs);

audiowrite('assignment9_bpf_source1.wav', sep1, fs);
audiowrite('assignment9_bpf_source2.wav', sep2, fs);

disp('Assignment 9 completed successfully.');
disp('Please take screenshots of the 8 figures manually.');
disp('Audio files saved.');
