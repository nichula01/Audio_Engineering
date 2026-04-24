% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 08
% Subspace Filtering
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script performs PCA-based subspace filtering on a noisy
% guitar signal using two sliding-window sizes:
%   1) window size = 2
%   2) window size = 20
%
% The idea is to keep only the dominant pattern (1st principal
% component) and reconstruct an approximation of the signal.
%
% NOTE:
% - Figures are shown on screen
% - Please take screenshots manually in MATLAB Online
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  Locate uploaded guitar file automatically
% --------------------------------------------------------------
mp3Files = dir('*guitar*.mp3');
wavFiles = dir('*guitar*.wav');

if ~isempty(mp3Files)
    audioFileName = mp3Files(1).name;
elseif ~isempty(wavFiles)
    audioFileName = wavFiles(1).name;
else
    error('No guitar audio file found. Please upload a file containing "guitar" in its name.');
end

disp(['Using audio file: ', audioFileName]);

%% -------------------------------------------------------------
%  Load and prepare audio
% --------------------------------------------------------------
[fileData, Fs] = audioread(audioFileName);

if size(fileData,2) > 1
    x = mean(fileData,2);
else
    x = fileData;
end

x = x(:);
x = x / max(abs(x) + eps);

% Use a short segment for faster processing
start_sec = 0;
dur_sec   = 3;

start_sample = floor(start_sec * Fs) + 1;
end_sample   = min(length(x), floor((start_sec + dur_sec) * Fs));

x = x(start_sample:end_sample);
x = x / max(abs(x) + eps);

N = length(x);
t = (0:N-1)/Fs;

%% -------------------------------------------------------------
%  Add white Gaussian noise manually (10 dB SNR)
% --------------------------------------------------------------
SNR_dB = 10;

signalPower = mean(x.^2);
noisePower  = signalPower / (10^(SNR_dB/10));

noise = sqrt(noisePower) * randn(size(x));
x_noisy = x + noise;
x_noisy = x_noisy / max(abs(x_noisy) + eps);

%% -------------------------------------------------------------
%  Apply PCA-based subspace filtering for two window sizes
% --------------------------------------------------------------
window_size_1 = 2;
window_size_2 = 20;

[X1, X1_centered, mean1] = create_sliding_matrix(x_noisy, window_size_1);
[X2, X2_centered, mean2] = create_sliding_matrix(x_noisy, window_size_2);

% PCA for window size 2
R1 = cov(X1_centered);
[V1, D1] = eig(R1);
[evals1, idx1] = sort(diag(D1), 'descend');
V1 = V1(:, idx1);

Y1 = X1_centered * V1;
Y1_pc1 = Y1(:,1);

X1_approx_centered = Y1_pc1 * V1(:,1)';
X1_approx = X1_approx_centered + mean1;

x1_filtered = reconstruct_from_windows(X1_approx);
x1_filtered = x1_filtered / max(abs(x1_filtered) + eps);

variance_explained_1 = evals1 / sum(evals1 + eps);

% PCA for window size 20
R2 = cov(X2_centered);
[V2, D2] = eig(R2);
[evals2, idx2] = sort(diag(D2), 'descend');
V2 = V2(:, idx2);

Y2 = X2_centered * V2;
Y2_pc1 = Y2(:,1);

X2_approx_centered = Y2_pc1 * V2(:,1)';
X2_approx = X2_approx_centered + mean2;

x2_filtered = reconstruct_from_windows(X2_approx);
x2_filtered = x2_filtered / max(abs(x2_filtered) + eps);

variance_explained_2 = evals2 / sum(evals2 + eps);

%% -------------------------------------------------------------
%  Prepare reduced vectors for cleaner plotting
% --------------------------------------------------------------
step_t = max(1, floor(length(x_noisy) / 2000));
idx_t  = 1:step_t:length(x_noisy);

t_plot       = t(idx_t);
x_noisy_plot = x_noisy(idx_t);
x1_plot      = x1_filtered(idx_t);
x2_plot      = x2_filtered(idx_t);

% Frequency-domain preparation (first half only)
Nf = length(x_noisy);
halfN = floor(Nf/2);

f_axis = (0:halfN-1) * (Fs/Nf);

X_noisy_fft = abs(fft(x_noisy));
X1_fft      = abs(fft(x1_filtered));
X2_fft      = abs(fft(x2_filtered));

X_noisy_half = X_noisy_fft(1:halfN);
X1_half      = X1_fft(1:halfN);
X2_half      = X2_fft(1:halfN);

step_f = max(1, floor(halfN / 2500));
idx_f  = 1:step_f:halfN;

f_plot       = f_axis(idx_f);
X_noisy_plot = X_noisy_half(idx_f);
X1_spec_plot = X1_half(idx_f);
X2_spec_plot = X2_half(idx_f);

%% -------------------------------------------------------------
%  Figure 1: Phase Space Analysis
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1);
plot(X1(:,1), X1(:,2), '.', 'MarkerSize', 2);
xlabel('x(i)');
ylabel('x(i+1)');
title('Phase Plot (Window Size = 2)');
grid on;

subplot(1,2,2);
plot(X2(:,1), X2(:,2), '.', 'MarkerSize', 2);
xlabel('x(i)');
ylabel('x(i+1)');
title('Phase Plot (Window Size = 20)');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 2: Projection onto First Principal Component
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1);
plot(Y1_pc1, 'LineWidth', 1);
xlabel('Sample Index');
ylabel('Amplitude');
title('Projection onto 1st PC (Window Size = 2)');
grid on;

subplot(1,2,2);
plot(Y2_pc1, 'LineWidth', 1);
xlabel('Sample Index');
ylabel('Amplitude');
title('Projection onto 1st PC (Window Size = 20)');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 3: Reconstructed Phase Space
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1);
plot(X1_approx(:,1), X1_approx(:,2), '.', 'MarkerSize', 2);
xlabel('Reconstructed x(i)');
ylabel('Reconstructed x(i+1)');
title('Reconstructed Phase Plot (Window Size = 2)');
grid on;

subplot(1,2,2);
plot(X2_approx(:,1), X2_approx(:,2), '.', 'MarkerSize', 2);
xlabel('Reconstructed x(i)');
ylabel('Reconstructed x(i+1)');
title('Reconstructed Phase Plot (Window Size = 20)');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Figure 4: Time-Domain Comparison
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1);
plot(t_plot, x_noisy_plot, 'b'); hold on;
plot(t_plot, x1_plot, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy vs Reconstructed (Window Size = 2)');
legend('Noisy Signal', 'Reconstructed');
grid on;
xlim([0 min(0.12, t(end))]);

subplot(1,2,2);
plot(t_plot, x_noisy_plot, 'b'); hold on;
plot(t_plot, x2_plot, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Noisy vs Reconstructed (Window Size = 20)');
legend('Noisy Signal', 'Reconstructed');
grid on;
xlim([0 min(0.12, t(end))]);

drawnow;

%% -------------------------------------------------------------
%  Figure 5: Frequency-Domain Comparison
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1100 450]);

subplot(1,2,1);
plot(f_plot, X_noisy_plot, 'b'); hold on;
plot(f_plot, X1_spec_plot, 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum Comparison (Window Size = 2)');
legend('Noisy Signal', 'Reconstructed');
grid on;
xlim([0 Fs/2]);

subplot(1,2,2);
plot(f_plot, X_noisy_plot, 'b'); hold on;
plot(f_plot, X2_spec_plot, 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Spectrum Comparison (Window Size = 20)');
legend('Noisy Signal', 'Reconstructed');
grid on;
xlim([0 Fs/2]);

drawnow;

%% -------------------------------------------------------------
%  Display variance explained
% --------------------------------------------------------------
disp('Variance explained by 1st principal component:');
disp(['Window size 2  : ', num2str(variance_explained_1(1)*100, '%.2f'), '%']);
disp(['Window size 20 : ', num2str(variance_explained_2(1)*100, '%.2f'), '%']);

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment8_original_segment.wav', x, Fs);
audiowrite('assignment8_noisy_signal.wav', x_noisy, Fs);
audiowrite('assignment8_reconstructed_ws2.wav', x1_filtered, Fs);
audiowrite('assignment8_reconstructed_ws20.wav', x2_filtered, Fs);

disp('Assignment 8 completed successfully.');
disp('Please take screenshots of the 5 figures manually.');
disp('Audio files saved:');
disp('  assignment8_original_segment.wav');
disp('  assignment8_noisy_signal.wav');
disp('  assignment8_reconstructed_ws2.wav');
disp('  assignment8_reconstructed_ws20.wav');

%% =============================================================
%  Local functions
% =============================================================
function [X, X_centered, meanVec] = create_sliding_matrix(signal, windowSize)
    signal = signal(:);
    N = length(signal);
    numRows = N - windowSize + 1;

    X = zeros(numRows, windowSize);
    for k = 1:windowSize
        X(:,k) = signal(k:k+numRows-1);
    end

    meanVec = mean(X, 1);
    X_centered = X - meanVec;
end

function x_rec = reconstruct_from_windows(Xwin)
    [numRows, windowSize] = size(Xwin);
    N = numRows + windowSize - 1;

    x_rec = zeros(N,1);
    counts = zeros(N,1);

    for i = 1:numRows
        idx = i:(i+windowSize-1);
        x_rec(idx) = x_rec(idx) + Xwin(i,:)';
        counts(idx) = counts(idx) + 1;
    end

    x_rec = x_rec ./ (counts + eps);
end
