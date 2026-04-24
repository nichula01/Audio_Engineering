% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 13
% Wavelet Transform (WT) / Scalogram
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Topic:
% Compare Sitar and Flute clips using:
%   1) Spectrograms
%   2) Continuous Wavelet Transform (CWT Scalogram)
%   3) Haar DWT
%   4) Daubechies db4 DWT
%
% MATLAB Online friendly:
% - figures shown on screen
% - no automatic image export
% - take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  File names
% --------------------------------------------------------------
sitarFile = 'Sitar_track.mp3';
fluteFile = 'Flute_track.mp3';

if ~isfile(sitarFile)
    error('File not found: Sitar_track.mp3');
end

if ~isfile(fluteFile)
    error('File not found: Flute_track.mp3');
end

disp(['Using file: ', sitarFile]);
disp(['Using file: ', fluteFile]);

%% -------------------------------------------------------------
%  Load audio
% --------------------------------------------------------------
[sitar, fsS] = audioread(sitarFile);
[flute, fsF] = audioread(fluteFile);

% Convert stereo to mono
if size(sitar,2) > 1
    sitar = mean(sitar,2);
end

if size(flute,2) > 1
    flute = mean(flute,2);
end

%% -------------------------------------------------------------
%  Trim to first 15 seconds (or shorter if clip is short)
% --------------------------------------------------------------
maxSamplesS = min(length(sitar), 15*fsS);
maxSamplesF = min(length(flute), 15*fsF);

sitar = sitar(1:maxSamplesS);
flute = flute(1:maxSamplesF);

% Normalize
sitar = sitar / max(abs(sitar) + eps);
flute = flute / max(abs(flute) + eps);

tS = (0:length(sitar)-1)/fsS;
tF = (0:length(flute)-1)/fsF;

%% -------------------------------------------------------------
%  Figure 1: Spectrogram comparison
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 1000 700]);

subplot(2,1,1);
spectrogram(sitar,256,200,256,fsS,'yaxis');
title('Spectrogram - Sitar');
colorbar;

subplot(2,1,2);
spectrogram(flute,256,200,256,fsF,'yaxis');
title('Spectrogram - Flute');
colorbar;

drawnow;

%% -------------------------------------------------------------
%  Figure 2: 2D CWT scalogram comparison
% --------------------------------------------------------------
[cfsS, frqS] = cwt(sitar, fsS, 'morse');
[cfsF, frqF] = cwt(flute, fsF, 'morse');

figure('Color','w','Position',[100 100 1000 700]);

subplot(2,1,1);
imagesc(tS, frqS, abs(cfsS));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('CWT Scalogram - Sitar');
colorbar;

subplot(2,1,2);
imagesc(tF, frqF, abs(cfsF));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('CWT Scalogram - Flute');
colorbar;

drawnow;

%% -------------------------------------------------------------
%  DWT parameters
% --------------------------------------------------------------
levels = 2;

%% -------------------------------------------------------------
%  Haar DWT
% --------------------------------------------------------------
[cS_haar, lS_haar] = wavedec(sitar, levels, 'haar');
[cF_haar, lF_haar] = wavedec(flute, levels, 'haar');

for lvl = 1:levels
    approxS = appcoef(cS_haar, lS_haar, 'haar', lvl);
    approxF = appcoef(cF_haar, lF_haar, 'haar', lvl);

    figure('Color','w','Position',[100 100 1000 500]);

    subplot(2,1,1);
    plot(approxS, 'LineWidth', 1);
    title(sprintf('Sitar - Haar DWT Level %d', lvl));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;

    subplot(2,1,2);
    plot(approxF, 'LineWidth', 1);
    title(sprintf('Flute - Haar DWT Level %d', lvl));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;

    drawnow;
end

%% -------------------------------------------------------------
%  Daubechies db4 DWT
% --------------------------------------------------------------
[cS_db4, lS_db4] = wavedec(sitar, levels, 'db4');
[cF_db4, lF_db4] = wavedec(flute, levels, 'db4');

for lvl = 1:levels
    approxS = appcoef(cS_db4, lS_db4, 'db4', lvl);
    approxF = appcoef(cF_db4, lF_db4, 'db4', lvl);

    figure('Color','w','Position',[100 100 1000 500]);

    subplot(2,1,1);
    plot(approxS, 'LineWidth', 1);
    title(sprintf('Sitar - db4 DWT Level %d', lvl));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;

    subplot(2,1,2);
    plot(approxF, 'LineWidth', 1);
    title(sprintf('Flute - db4 DWT Level %d', lvl));
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;

    drawnow;
end

%% -------------------------------------------------------------
%  Save audio clips
% --------------------------------------------------------------
audiowrite('assignment13_sitar_clip.wav', sitar, fsS);
audiowrite('assignment13_flute_clip.wav', flute, fsF);

disp('Assignment 13 completed successfully.');
disp('Please take screenshots manually.');
disp('Suggested screenshot names:');
disp('  assignment13_spectrogram_comparison.png');
disp('  assignment13_cwt_scalogram.png');
disp('  assignment13_haar_level1.png');
disp('  assignment13_haar_level2.png');
disp('  assignment13_db4_level1.png');
disp('  assignment13_db4_level2.png');
disp('Audio files saved:');
disp('  assignment13_sitar_clip.wav');
disp('  assignment13_flute_clip.wav');
