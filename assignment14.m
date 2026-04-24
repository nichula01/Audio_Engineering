% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 14
% Mel-Frequency Cepstral Coefficients (MFCC)
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Description:
% Compare MFCC features of two different sounds using:
%   1) Pre-emphasis PSD
%   2) Log-Mel spectrogram
%   3) MFCC heatmap
%   4) Mean MFCC comparison
%   5) 3D MFCC scatter
%   6) PCA projection
%
% Toolbox-friendly note:
% - Does NOT require MATLAB mfcc() function
% - Uses manual MFCC extraction
% - Figures shown on screen
% - Please take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  File names
% --------------------------------------------------------------
file1 = 'Sound1_track.mp3';
file2 = 'Sound2_track.mp3';

if ~isfile(file1)
    error('File not found: Sound1_track.mp3');
end
if ~isfile(file2)
    error('File not found: Sound2_track.mp3');
end

disp(['Using file 1: ', file1]);
disp(['Using file 2: ', file2]);

%% -------------------------------------------------------------
%  Load audio
% --------------------------------------------------------------
[x1, Fs1] = audioread(file1);
[x2, Fs2] = audioread(file2);

if size(x1,2) > 1
    x1 = mean(x1,2);
end
if size(x2,2) > 1
    x2 = mean(x2,2);
end

%% -------------------------------------------------------------
%  Resample if needed
% --------------------------------------------------------------
Fs = 16000;
if Fs1 ~= Fs
    x1 = resample(x1, Fs, Fs1);
end
if Fs2 ~= Fs
    x2 = resample(x2, Fs, Fs2);
end

%% -------------------------------------------------------------
%  Trim to first 20 seconds for speed and stability
% --------------------------------------------------------------
x1 = x1(1:min(length(x1), 20*Fs));
x2 = x2(1:min(length(x2), 20*Fs));

x1 = x1 / max(abs(x1) + eps);
x2 = x2 / max(abs(x2) + eps);

%% -------------------------------------------------------------
%  Pre-emphasis
% --------------------------------------------------------------
alpha = 0.95;
x1 = filter([1 -alpha], 1, x1);
x2 = filter([1 -alpha], 1, x2);

%% -------------------------------------------------------------
%  Pre-emphasis PSD (Welch)
% --------------------------------------------------------------
[P1, f1_psd] = pwelch(x1, [], [], [], Fs);
[P2, f2_psd] = pwelch(x2, [], [], [], Fs);

figure('Color','w','Position',[100 100 900 450]);
plot(f1_psd, 10*log10(P1 + eps), 'b', 'LineWidth', 1.2); hold on;
plot(f2_psd, 10*log10(P2 + eps), 'r', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Welch PSD After Pre-emphasis');
legend('Sound 1','Sound 2');
grid on;
drawnow;

%% -------------------------------------------------------------
%  Framing parameters
% --------------------------------------------------------------
frameLen   = round(0.025 * Fs);   % 25 ms
frameShift = round(0.010 * Fs);   % 10 ms
NFFT       = 512;
numCoeffs  = 13;
numFilters = 26;

win = hamming(frameLen, 'periodic');

%% -------------------------------------------------------------
%  Manual MFCC extraction for both signals
% --------------------------------------------------------------
[coeffs1, logMel1] = manual_mfcc(x1, Fs, frameLen, frameShift, win, NFFT, numFilters, numCoeffs);
[coeffs2, logMel2] = manual_mfcc(x2, Fs, frameLen, frameShift, win, NFFT, numFilters, numCoeffs);

%% -------------------------------------------------------------
%  Figure 2: Log-Mel spectrogram
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 500]);

subplot(1,2,1);
imagesc(logMel1');
axis xy;
title('Sound 1 - Log Mel Spectrogram');
xlabel('Frame Index');
ylabel('Mel Band');
colorbar;

subplot(1,2,2);
imagesc(logMel2');
axis xy;
title('Sound 2 - Log Mel Spectrogram');
xlabel('Frame Index');
ylabel('Mel Band');
colorbar;

drawnow;

%% -------------------------------------------------------------
%  Figure 3: MFCC heatmap
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 500]);

subplot(1,2,1);
imagesc(coeffs1');
axis xy;
title('Sound 1 - MFCC (C0 to C12)');
xlabel('Frame Index');
ylabel('Coefficient Index');
colorbar;

subplot(1,2,2);
imagesc(coeffs2');
axis xy;
title('Sound 2 - MFCC (C0 to C12)');
xlabel('Frame Index');
ylabel('Coefficient Index');
colorbar;

drawnow;

%% -------------------------------------------------------------
%  Remove C0
% --------------------------------------------------------------
coeffs1_noc0 = coeffs1(:,2:end);
coeffs2_noc0 = coeffs2(:,2:end);

%% -------------------------------------------------------------
%  Mean MFCC comparison
% --------------------------------------------------------------
mean1 = mean(coeffs1_noc0, 1);
mean2 = mean(coeffs2_noc0, 1);

figure('Color','w','Position',[100 100 900 450]);
plot(mean1, 'b-o', 'LineWidth', 1.2); hold on;
plot(mean2, 'r-o', 'LineWidth', 1.2);
xlabel('MFCC Coefficient Index');
ylabel('Mean Value');
title('Mean MFCC Comparison (Without C0)');
legend('Sound 1','Sound 2');
grid on;
drawnow;

%% -------------------------------------------------------------
%  Similarity measures
% --------------------------------------------------------------
euclidDist = norm(mean1 - mean2);
cosSim = dot(mean1, mean2) / (norm(mean1)*norm(mean2) + eps);

% Simple DTW on MFCC sequences
dtwDist = simple_dtw(coeffs1_noc0', coeffs2_noc0');

fprintf('\n---- MFCC Comparison Results ----\n');
fprintf('Euclidean Distance : %.4f\n', euclidDist);
fprintf('Cosine Similarity  : %.4f\n', cosSim);
fprintf('DTW Distance       : %.4f\n', dtwDist);

%% -------------------------------------------------------------
%  Figure 5: 3D MFCC scatter (C1,C2,C3)
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 500]);
scatter3(coeffs1_noc0(:,1), coeffs1_noc0(:,2), coeffs1_noc0(:,3), ...
    10, 'b', 'filled'); hold on;
scatter3(coeffs2_noc0(:,1), coeffs2_noc0(:,2), coeffs2_noc0(:,3), ...
    10, 'r', 'filled');
xlabel('C1');
ylabel('C2');
zlabel('C3');
title('3D MFCC Scatter');
legend('Sound 1','Sound 2');
grid on;
view(40,30);
drawnow;

%% -------------------------------------------------------------
%  Figure 6: PCA projection
% --------------------------------------------------------------
allFeatures = [coeffs1_noc0; coeffs2_noc0];
labels = [ones(size(coeffs1_noc0,1),1); 2*ones(size(coeffs2_noc0,1),1)];

% Manual PCA using SVD
X = allFeatures - mean(allFeatures,1);
[U,S,~] = svd(X, 'econ');
score = U*S;

varvals = diag(S).^2;
explained = 100 * varvals / sum(varvals + eps);
expl3 = sum(explained(1:min(3,end)));

figure('Color','w','Position',[100 100 900 500]);
scatter3(score(labels==1,1), score(labels==1,2), score(labels==1,3), ...
    10, 'b', 'filled'); hold on;
scatter3(score(labels==2,1), score(labels==2,2), score(labels==2,3), ...
    10, 'r', 'filled');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title(sprintf('PCA Projection (Variance Explained: %.2f%%)', expl3));
legend('Sound 1','Sound 2');
grid on;
view(40,30);
drawnow;

%% -------------------------------------------------------------
%  Save processed audio
% --------------------------------------------------------------
audiowrite('assignment14_sound1_preemphasized.wav', x1 / max(abs(x1)+eps), Fs);
audiowrite('assignment14_sound2_preemphasized.wav', x2 / max(abs(x2)+eps), Fs);

disp('Assignment 14 completed successfully.');
disp('Please take screenshots manually.');
disp('Suggested screenshot names:');
disp('  assignment14_preemphasis_psd.png');
disp('  assignment14_logmel.png');
disp('  assignment14_mfcc_heatmap.png');
disp('  assignment14_mean_mfcc.png');
disp('  assignment14_mfcc_3d.png');
disp('  assignment14_pca_projection.png');
disp('Audio files saved:');
disp('  assignment14_sound1_preemphasized.wav');
disp('  assignment14_sound2_preemphasized.wav');

%% =============================================================
%  Local functions
% =============================================================
function [coeffs, logMel] = manual_mfcc(x, Fs, frameLen, frameShift, win, NFFT, numFilters, numCoeffs)
    x = x(:);

    % Frame using buffer-style manual indexing
    numFrames = 1 + floor((length(x) - frameLen) / frameShift);
    frames = zeros(numFrames, frameLen);

    for i = 1:numFrames
        idx = (1:frameLen) + (i-1)*frameShift;
        frames(i,:) = x(idx);
    end

    % Windowing
    frames = frames .* win';

    % Power spectrum
    spec = abs(fft(frames, NFFT, 2)).^2;
    spec = spec(:, 1:NFFT/2+1);

    % Mel filterbank
    fb = mel_filterbank(numFilters, NFFT, Fs);

    % Mel energies
    melEnergy = spec * fb';
    logMel = log(melEnergy + 1e-6);

    % DCT
    coeffs_all = dct(logMel, [], 2);
    coeffs = coeffs_all(:, 1:numCoeffs);
end

function fb = mel_filterbank(M, NFFT, Fs)
    fmin = 0;
    fmax = Fs/2;

    melmin = 2595*log10(1 + fmin/700);
    melmax = 2595*log10(1 + fmax/700);

    melpts = linspace(melmin, melmax, M+2);
    hzpts = 700*(10.^(melpts/2595) - 1);

    bins = floor((NFFT+1)*hzpts/Fs);

    fb = zeros(M, NFFT/2+1);

    for m = 2:M+1
        f_left  = bins(m-1);
        f_center= bins(m);
        f_right = bins(m+1);

        for k = f_left:f_center
            if k >= 0 && k <= NFFT/2
                fb(m-1, k+1) = (k - f_left) / max(f_center - f_left, 1);
            end
        end

        for k = f_center:f_right
            if k >= 0 && k <= NFFT/2
                fb(m-1, k+1) = (f_right - k) / max(f_right - f_center, 1);
            end
        end
    end
end

function dist = simple_dtw(A, B)
    % A and B are feature_dim x time_frames
    na = size(A,2);
    nb = size(B,2);

    D = inf(na+1, nb+1);
    D(1,1) = 0;

    for i = 2:na+1
        for j = 2:nb+1
            cost = norm(A(:,i-1) - B(:,j-1));
            D(i,j) = cost + min([D(i-1,j), D(i,j-1), D(i-1,j-1)]);
        end
    end

    dist = D(end,end);
end
