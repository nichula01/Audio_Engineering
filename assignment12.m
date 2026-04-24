% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 12
% Short Time Fourier Transform (STFT) / Spectrogram
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Description:
% This script mixes two instrument signals (guitar and sansula),
% analyzes the mixture using FT and STFT, performs simple
% frequency-mask separation, and reconstructs the separated audio.
%
% MATLAB Online friendly:
% - figures shown on screen
% - no auto image export
% - please take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  Locate audio files automatically
% --------------------------------------------------------------
guitar_mp3 = dir('*guitar*.mp3');
guitar_wav = dir('*guitar*.wav');

sansula_mp3 = dir('*sansula*.mp3');
sansula_wav = dir('*sansula*.wav');

if ~isempty(guitar_mp3)
    guitarFile = guitar_mp3(1).name;
elseif ~isempty(guitar_wav)
    guitarFile = guitar_wav(1).name;
else
    error('No guitar audio file found.');
end

if ~isempty(sansula_mp3)
    sansulaFile = sansula_mp3(1).name;
elseif ~isempty(sansula_wav)
    sansulaFile = sansula_wav(1).name;
else
    error('No sansula audio file found.');
end

disp(['Using guitar file  : ', guitarFile]);
disp(['Using sansula file : ', sansulaFile]);

%% -------------------------------------------------------------
%  Load audio
% --------------------------------------------------------------
[song1_raw, fs1] = audioread(guitarFile);
[song2_raw, fs2] = audioread(sansulaFile);

if size(song1_raw,2) > 1, song1_raw = mean(song1_raw,2); end
if size(song2_raw,2) > 1, song2_raw = mean(song2_raw,2); end

%% -------------------------------------------------------------
%  Use a common sample rate
% --------------------------------------------------------------
fs = 8000;

if fs1 ~= fs
    song1_raw = resample(song1_raw, fs, fs1);
end

if fs2 ~= fs
    song2_raw = resample(song2_raw, fs, fs2);
end

%% -------------------------------------------------------------
%  Trim to same length and shorten for stable MATLAB Online use
% --------------------------------------------------------------
len = min(length(song1_raw), length(song2_raw));
song1 = song1_raw(1:len);
song2 = song2_raw(1:len);

maxDur = 4;  % seconds
maxLen = min(length(song1), maxDur*fs);

song1 = song1(1:maxLen);
song2 = song2(1:maxLen);

song1 = song1 / max(abs(song1) + eps);
song2 = song2 / max(abs(song2) + eps);

t = (0:length(song1)-1)/fs;

%% -------------------------------------------------------------
%  Mixture
% --------------------------------------------------------------
mix = song1 + song2;
mix = mix / max(abs(mix) + eps);

%% -------------------------------------------------------------
%  Fourier Transform of mixture
% --------------------------------------------------------------
N = length(mix);
Mix_FT = abs(fft(mix));
halfN = floor(N/2);
f = (0:halfN-1)*(fs/N);

%% -------------------------------------------------------------
%  STFT parameters
% --------------------------------------------------------------
window  = 256;
overlap = 192;
nfft    = 512;
win     = hamming(window, 'periodic');

[S_mix,  F,  T ] = spectrogram(mix,   win, overlap, nfft, fs);
[S1_orig,F1, T1] = spectrogram(song1, win, overlap, nfft, fs);
[S2_orig,F2, T2] = spectrogram(song2, win, overlap, nfft, fs);

%% -------------------------------------------------------------
%  Frequency-mask separation
% --------------------------------------------------------------
% Tune this if needed after looking at the mixture spectrum
cutoff = 1000;   % Hz

mask1 = F < cutoff;      % lower-frequency region
mask2 = F >= cutoff;     % higher-frequency region

S1_ext = S_mix .* mask1;
S2_ext = S_mix .* mask2;

%% -------------------------------------------------------------
%  Reconstruct time-domain separated signals
% --------------------------------------------------------------
% Use inverse STFT if available
canUseIstft = exist('istft', 'file') == 2;

if canUseIstft
    song1_rec = istft(S1_ext, fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
    song2_rec = istft(S2_ext, fs, 'Window', win, 'OverlapLength', overlap, 'FFTLength', nfft);
else
    % Fallback: keep outputs zero-sized placeholders if istft is unavailable
    song1_rec = zeros(size(mix));
    song2_rec = zeros(size(mix));
end

% Match lengths safely
Lrec1 = min(length(song1_rec), length(mix));
Lrec2 = min(length(song2_rec), length(mix));

song1_rec = song1_rec(1:Lrec1);
song2_rec = song2_rec(1:Lrec2);

if max(abs(song1_rec)) > 0
    song1_rec = song1_rec / max(abs(song1_rec) + eps);
end
if max(abs(song2_rec)) > 0
    song2_rec = song2_rec / max(abs(song2_rec) + eps);
end

%% -------------------------------------------------------------
%  Reduced vectors for plotting
% --------------------------------------------------------------
step_t = max(1, floor(length(mix)/2000));
idx_t  = 1:step_t:length(mix);

t_plot   = t(idx_t);
mix_plot = mix(idx_t);

%% -------------------------------------------------------------
%  Figure 1: Fourier Transform of Mixture
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 420]);
plot(f, Mix_FT(1:halfN), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Fourier Transform of Mixture Signal');
grid on;
xlim([0 2000]);
drawnow;

%% -------------------------------------------------------------
%  Figure 2: STFT Spectrogram of Mixture
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 900 500]);
imagesc(T, F, 20*log10(abs(S_mix) + 1e-6));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('STFT Spectrogram of Mixture');
colorbar;
ylim([0 2500]);
drawnow;

%% -------------------------------------------------------------
%  Figure 3: Original vs Extracted Song 1
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 950 650]);

subplot(2,1,1);
imagesc(T1, F1, 20*log10(abs(S1_orig) + 1e-6));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original Guitar Spectrogram');
colorbar;
ylim([0 2500]);

subplot(2,1,2);
imagesc(T, F, 20*log10(abs(S1_ext) + 1e-6));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Extracted Guitar from Mixture');
colorbar;
ylim([0 2500]);

drawnow;

%% -------------------------------------------------------------
%  Figure 4: Original vs Extracted Song 2
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 950 650]);

subplot(2,1,1);
imagesc(T2, F2, 20*log10(abs(S2_orig) + 1e-6));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Original Sansula Spectrogram');
colorbar;
ylim([0 2500]);

subplot(2,1,2);
imagesc(T, F, 20*log10(abs(S2_ext) + 1e-6));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Extracted Sansula from Mixture');
colorbar;
ylim([0 2500]);

drawnow;

%% -------------------------------------------------------------
%  Figure 5: Reconstructed separated signals (time domain)
% --------------------------------------------------------------
if canUseIstft
    figure('Color','w','Position',[100 100 950 500]);

    subplot(2,1,1);
    plot((0:length(song1_rec)-1)/fs, song1_rec, 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Recovered Guitar (Time Domain)');
    grid on;

    subplot(2,1,2);
    plot((0:length(song2_rec)-1)/fs, song2_rec, 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Recovered Sansula (Time Domain)');
    grid on;

    drawnow;
end

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment12_guitar_original.wav', song1, fs);
audiowrite('assignment12_sansula_original.wav', song2, fs);
audiowrite('assignment12_mixture.wav', mix, fs);

if canUseIstft
    audiowrite('assignment12_guitar_recovered.wav', song1_rec, fs);
    audiowrite('assignment12_sansula_recovered.wav', song2_rec, fs);
end

disp('Assignment 12 completed successfully.');
disp('Please take screenshots manually.');
disp('Suggested screenshot names:');
disp('  assignment12_mixture_ft.png');
disp('  assignment12_mixture_stft.png');
disp('  assignment12_guitar_extraction.png');
disp('  assignment12_sansula_extraction.png');
if canUseIstft
    disp('  assignment12_recovered_time_domain.png');
end

disp('Audio files saved:');
disp('  assignment12_guitar_original.wav');
disp('  assignment12_sansula_original.wav');
disp('  assignment12_mixture.wav');
if canUseIstft
    disp('  assignment12_guitar_recovered.wav');
    disp('  assignment12_sansula_recovered.wav');
end
