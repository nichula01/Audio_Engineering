% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 15
% Audio Signal Compression and Reconstruction using DCT
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Description:
% This script compresses an audio signal using Discrete Cosine
% Transform (DCT) coefficients and reconstructs it using
% 5%, 10%, and 20% retained coefficients.
%
% Outputs:
%   - Original signal and DCT coefficient magnitude
%   - Time-domain reconstruction comparison
%   - Frequency-domain comparison
%   - Cumulative DCT energy plot
%   - Metrics: MSE, SNR, Compression Ratio
%   - Reconstructed WAV files
%
% MATLAB Online friendly:
% - figures shown on screen
% - no automatic image export
% - please take screenshots manually
% ===============================================================

clc;
clear;
close all;

set(groot, 'defaultFigureVisible', 'on');

%% -------------------------------------------------------------
%  Locate audio file
% --------------------------------------------------------------
if isfile('Piano_Track.mp3')
    audioFileName = 'Piano_Track.mp3';
elseif isfile('Piano_Track.wav')
    audioFileName = 'Piano_Track.wav';
else
    mp3Files = dir('*.mp3');
    wavFiles = dir('*.wav');

    if ~isempty(mp3Files)
        audioFileName = mp3Files(1).name;
    elseif ~isempty(wavFiles)
        audioFileName = wavFiles(1).name;
    else
        error('No audio file found. Upload Piano_Track.mp3 or another audio file.');
    end
end

disp(['Using audio file: ', audioFileName]);

%% -------------------------------------------------------------
%  Load audio
% --------------------------------------------------------------
[fileData, Fs] = audioread(audioFileName);

if size(fileData,2) == 2
    fileData = mean(fileData,2);
end

x = fileData(:);
x = x / max(abs(x) + eps);

% Use up to 12 seconds for stable MATLAB Online operation
maxLen = min(length(x), 12*Fs);
x = x(1:maxLen);

N = length(x);
t = (0:N-1)/Fs;

fprintf('Original Signal Length: %d samples\n', N);

%% -------------------------------------------------------------
%  Apply DCT
% --------------------------------------------------------------
X = dct(x);

%% -------------------------------------------------------------
%  Reduced vectors for plotting
% --------------------------------------------------------------
step_t = max(1, floor(N / 2500));
idx_t  = 1:step_t:N;

t_plot = t(idx_t);
x_plot = x(idx_t);

% DCT magnitude for display
Xmag = abs(X);
step_c = max(1, floor(N / 2500));
idx_c  = 1:step_c:N;

coeff_plot_idx = idx_c;
Xmag_plot = Xmag(idx_c);

%% -------------------------------------------------------------
%  Figure 1: Original signal and DCT coefficients
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 950 650]);

subplot(2,1,1);
plot(t_plot, x_plot, 'b', 'LineWidth', 1);
title('Original Audio Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal');
grid on;

subplot(2,1,2);
plot(coeff_plot_idx, Xmag_plot, 'r', 'LineWidth', 1);
title('Magnitude of DCT Coefficients');
xlabel('Coefficient Index');
ylabel('Magnitude');
legend('DCT Magnitude');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Compression percentages
% --------------------------------------------------------------
retainPercents = [5 10 20];

reconstructedSignals = cell(length(retainPercents),1);
metrics = zeros(length(retainPercents), 3); 
% columns: [MSE, SNR, CompressionRatio]

%% -------------------------------------------------------------
%  Reconstruct using top-magnitude DCT coefficients
% --------------------------------------------------------------
for i = 1:length(retainPercents)
    retainPercent = retainPercents(i);

    numKeep = round((retainPercent/100) * N);

    % Keep largest-magnitude coefficients
    [~, idxSorted] = sort(abs(X), 'descend');
    keepIdx = idxSorted(1:numKeep);

    Xc = zeros(size(X));
    Xc(keepIdx) = X(keepIdx);

    % Reconstruct
    x_rec = idct(Xc);
    x_rec = x_rec / max(abs(x_rec) + eps);

    reconstructedSignals{i} = x_rec;

    % Metrics
    mseVal = mean((x - x_rec).^2);
    snrVal = 10*log10(sum(x.^2) / (sum((x - x_rec).^2) + eps));
    compressionRatio = N / numKeep;

    metrics(i,:) = [mseVal, snrVal, compressionRatio];

    fprintf('\n==== %d%% Coefficients Retained ====\n', retainPercent);
    fprintf('Coefficients used: %d\n', numKeep);
    fprintf('MSE: %.6f\n', mseVal);
    fprintf('SNR: %.2f dB\n', snrVal);
    fprintf('Compression Ratio: %.2f : 1\n', compressionRatio);
end

%% -------------------------------------------------------------
%  Figure 2: Time-domain reconstruction comparison
% --------------------------------------------------------------
figure('Color','w','Position',[100 100 950 800]);

for i = 1:length(retainPercents)
    x_rec = reconstructedSignals{i};
    x_rec_plot = x_rec(idx_t);

    subplot(length(retainPercents),1,i);
    plot(t_plot, x_plot, 'b', 'LineWidth', 1); hold on;
    plot(t_plot, x_rec_plot, 'r', 'LineWidth', 1);
    title(sprintf('Time-Domain Reconstruction (%d%% Coefficients Retained)', retainPercents(i)));
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original', 'Reconstructed');
    grid on;
    xlim([0 min(0.12, t(end))]);
end

drawnow;

%% -------------------------------------------------------------
%  Figure 3: Frequency-domain comparison
% --------------------------------------------------------------
Xf = abs(fft(x));
halfN = floor(N/2);
f = (0:halfN-1)*(Fs/N);

Xf_half = Xf(1:halfN);
step_f = max(1, floor(halfN / 2500));
idx_f = 1:step_f:halfN;
f_plot = f(idx_f);
Xf_plot = Xf_half(idx_f);

figure('Color','w','Position',[100 100 950 800]);

for i = 1:length(retainPercents)
    x_rec = reconstructedSignals{i};
    Xrf = abs(fft(x_rec));
    Xrf_half = Xrf(1:halfN);
    Xrf_plot = Xrf_half(idx_f);

    subplot(length(retainPercents),1,i);
    plot(f_plot, Xf_plot, 'b', 'LineWidth', 1); hold on;
    plot(f_plot, Xrf_plot, 'r', 'LineWidth', 1);
    title(sprintf('Spectrum Comparison (%d%% Coefficients Retained)', retainPercents(i)));
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    legend('Original', 'Reconstructed');
    xlim([0 Fs/2]);
    grid on;
end

drawnow;

%% -------------------------------------------------------------
%  Figure 4: Cumulative DCT energy
% --------------------------------------------------------------
Xenergy = abs(X).^2;
XenergySorted = sort(Xenergy, 'descend');
cumEnergy = cumsum(XenergySorted) / sum(XenergySorted + eps);

figure('Color','w','Position',[100 100 900 450]);
plot((1:N)/N*100, cumEnergy*100, 'k', 'LineWidth', 1.5); hold on;
yline(5, '--r');
yline(10, '--b');
yline(20, '--g');
xlabel('Percentage of Coefficients Retained (%)');
ylabel('Cumulative Energy (%)');
title('Cumulative DCT Energy Retention');
grid on;

drawnow;

%% -------------------------------------------------------------
%  Save reconstructed audio
% --------------------------------------------------------------
audiowrite('assignment15_original_signal.wav', x, Fs);
audiowrite('assignment15_reconstructed_5.wav', reconstructedSignals{1}, Fs);
audiowrite('assignment15_reconstructed_10.wav', reconstructedSignals{2}, Fs);
audiowrite('assignment15_reconstructed_20.wav', reconstructedSignals{3}, Fs);

disp('Assignment 15 completed successfully.');
disp('Please take screenshots manually.');
disp('Suggested screenshot names:');
disp('  assignment15_original_and_dct.png');
disp('  assignment15_time_reconstruction.png');
disp('  assignment15_frequency_comparison.png');
disp('  assignment15_cumulative_energy.png');
disp('Audio files saved:');
disp('  assignment15_original_signal.wav');
disp('  assignment15_reconstructed_5.wav');
disp('  assignment15_reconstructed_10.wav');
disp('  assignment15_reconstructed_20.wav');
