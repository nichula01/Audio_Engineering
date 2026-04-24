% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 10
% Adaptive Wiener Filter
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script adds colored noise to a guitar signal and then
% applies two adaptive Wiener-style filters to recover the signal.
%
% Outputs:
%   - Original vs channel noise
%   - Original vs extracted (Case 1)
%   - Original vs extracted (Case 2)
%   - WAV files for listening
%
% NOTE:
% - MATLAB Online friendly
% - no automatic PNG export
% - take screenshots manually
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
    song = mean(fileData,2);
else
    song = fileData;
end

song = song(:);
song = song / max(abs(song) + eps);

% Use a shorter segment for speed
start_sec = 0;
dur_sec   = 3;

start_sample = floor(start_sec * Fs) + 1;
end_sample   = min(length(song), floor((start_sec + dur_sec) * Fs));

song = song(start_sample:end_sample);
song = song / max(abs(song) + eps);

N = length(song);
t = (0:N-1) / Fs;

%% -------------------------------------------------------------
%  Generate stronger white noise
% --------------------------------------------------------------
noise = 0.10 * randn(N,1);

%% -------------------------------------------------------------
%  Pass noise through a channel (color the noise)
% --------------------------------------------------------------
[bch, ach] = cheby1(4, 1, 0.4);     % 4th-order Chebyshev channel
noise_channel = filter(bch, ach, noise);
noise_channel = noise_channel / max(abs(noise_channel) + eps);

%% -------------------------------------------------------------
%  Observed signal = clean song + colored noise
% --------------------------------------------------------------
observed = song + 0.35 * noise_channel;
observed = observed / max(abs(observed) + eps);

%% -------------------------------------------------------------
%  Reduced vectors for plotting
% --------------------------------------------------------------
step_t = max(1, floor(N / 2000));
idx_t  = 1:step_t:N;

t_plot        = t(idx_t);
song_plot     = song(idx_t);
noise_plot    = noise_channel(idx_t);
observed_plot = observed(idx_t);

%% =============================================================
%  CASE 1: Weaker adaptive Wiener filter
% =============================================================
L1 = 8;                     % filter length
Q1 = 1e-3 * eye(L1);        % process noise covariance
R1 = 1e-1;                  % measurement noise covariance

w1 = zeros(L1,1);           % adaptive weights
P1 = eye(L1);               % error covariance
estimated1 = zeros(N,1);

for n = L1:N
    x_vec = observed(n:-1:n-L1+1);

    % prediction
    w_pred = w1;
    P_pred = P1 + Q1;

    % gain
    K = (P_pred * x_vec) / (x_vec' * P_pred * x_vec + R1);

    % estimate
    y_hat = w_pred' * x_vec;

    % error against clean reference
    e = song(n) - y_hat;

    % update
    w1 = w_pred + K * e;
    P1 = (eye(L1) - K * x_vec') * P_pred;

    estimated1(n) = w1' * x_vec;
end

estimated1 = estimated1 / max(abs(estimated1) + eps);

%% =============================================================
%  CASE 2: Stronger adaptive Wiener filter
% =============================================================
L2 = 20;                    % longer filter
Q2 = 1e-6 * eye(L2);        % smaller process noise
R2 = 1e-3;                  % smaller measurement noise

w2 = zeros(L2,1);
P2 = eye(L2);
estimated2 = zeros(N,1);

for n = L2:N
    x_vec = observed(n:-1:n-L2+1);

    % prediction
    w_pred = w2;
    P_pred = P2 + Q2;

    % gain
    K = (P_pred * x_vec) / (x_vec' * P_pred * x_vec + R2);

    % estimate
    y_hat = w_pred' * x_vec;

    % error against clean reference
    e = song(n) - y_hat;

    % update
    w2 = w_pred + K * e;
    P2 = (eye(L2) - K * x_vec') * P_pred;

    estimated2(n) = w2' * x_vec;
end

estimated2 = estimated2 / max(abs(estimated2) + eps);

%% -------------------------------------------------------------
%  Figure 1: Original signal and channel noise
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 450]);
plot(t_plot, song_plot, 'b', 'LineWidth', 1); hold on;
plot(t_plot, noise_plot, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Guitar Signal and Channel Noise');
legend('Original Signal', 'Noise after Channel');
grid on;
xlim([0 min(0.08, t(end))]);

drawnow;

%% -------------------------------------------------------------
%  Figure 2: Original vs Extracted (Case 1)
% --------------------------------------------------------------
estimated1_plot = estimated1(idx_t);

figure('Color','w','Position',[100 100 900 450]);
plot(t_plot, song_plot, 'b', 'LineWidth', 1); hold on;
plot(t_plot, estimated1_plot, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Extracted Signal (Case 1: Weaker Filter)');
legend('Original Signal', 'Extracted Signal - Case 1');
grid on;
xlim([0 min(0.08, t(end))]);

drawnow;

%% -------------------------------------------------------------
%  Figure 3: Original vs Extracted (Case 2)
% --------------------------------------------------------------
estimated2_plot = estimated2(idx_t);

figure('Color','w','Position',[100 100 900 450]);
plot(t_plot, song_plot, 'b', 'LineWidth', 1); hold on;
plot(t_plot, estimated2_plot, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Extracted Signal (Case 2: Stronger Filter)');
legend('Original Signal', 'Extracted Signal - Case 2');
grid on;
xlim([0 min(0.08, t(end))]);

drawnow;

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment10_original_segment.wav', song, Fs);
audiowrite('assignment10_observed_signal.wav', observed, Fs);
audiowrite('assignment10_extracted_case1.wav', estimated1, Fs);
audiowrite('assignment10_extracted_case2.wav', estimated2, Fs);

disp('Assignment 10 completed successfully.');
disp('Please take screenshots of the 3 figures manually.');
disp('Audio files saved:');
disp('  assignment10_original_segment.wav');
disp('  assignment10_observed_signal.wav');
disp('  assignment10_extracted_case1.wav');
disp('  assignment10_extracted_case2.wav');
