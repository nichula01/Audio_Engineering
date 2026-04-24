% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 05
% Filtering Music Signals
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script loads a guitar recording and applies:
%   1) Low-Pass Filter (LPF)
%   2) High-Pass Filter (HPF)
%   3) Band-Pass Filter (BPF)
%   4) Comb Filter
%
% For each filter, it generates:
%   - time-domain comparison
%   - filter frequency response
%   - frequency-domain comparison
%
% All figures are saved automatically as PNG files.
% ===============================================================

clc;
clear;
close all;

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

% Convert stereo to mono if needed
if size(fileData, 2) > 1
    x = mean(fileData, 2);
else
    x = fileData;
end

% Normalize full signal
x = x / max(abs(x) + eps);

% Select a short segment for analysis
start_sec = 0;
dur_sec   = 8;

start_sample = floor(start_sec * Fs) + 1;
end_sample   = min(length(x), floor((start_sec + dur_sec) * Fs));

x = x(start_sample:end_sample);

% Normalize trimmed signal
x = x / max(abs(x) + eps);

t = (0:length(x)-1) / Fs;

% FFT settings for later plots
Nfft = length(x);
fAxis = (0:Nfft-1) * (Fs / Nfft);

%% -------------------------------------------------------------
%  Common filter design settings
% --------------------------------------------------------------
Rp = 1;      % Passband ripple in dB
Rs = 50;     % Stopband attenuation in dB

%% =============================================================
%  1) LOW-PASS FILTER
% ==============================================================
Fp_lpf = 1800;      % Passband edge (Hz)
Fs_lpf = 2600;      % Stopband edge (Hz)

Wp_lpf = Fp_lpf / (Fs/2);
Ws_lpf = Fs_lpf / (Fs/2);

[n_lpf, Wn_lpf] = ellipord(Wp_lpf, Ws_lpf, Rp, Rs);
[b_lpf, a_lpf] = ellip(n_lpf, Rp, Rs, Wn_lpf, 'low');

y_lpf = filter(b_lpf, a_lpf, x);
y_lpf = y_lpf / max(abs(y_lpf) + eps);

% Frequency response
[H_lpf, f_lpf] = freqz(b_lpf, a_lpf, 2048, Fs);

% Spectra
X_mag    = abs(fft(x));
Y_lpfmag = abs(fft(y_lpf));

fig1 = figure;
subplot(3,1,1);
plot(t, x, 'b'); hold on;
plot(t, y_lpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Low-Pass Filter: Original vs Filtered (Time Domain)');
legend('Original','Filtered');
grid on;

subplot(3,1,2);
plot(f_lpf, abs(H_lpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Low-Pass Filter Frequency Response');
grid on;

subplot(3,1,3);
plot(fAxis, X_mag, 'b'); hold on;
plot(fAxis, Y_lpfmag, 'r');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Low-Pass Filter: Original vs Filtered (Frequency Domain)');
legend('Original','Filtered');
grid on;

saveas(fig1, 'assignment5_lpf_analysis.png');

%% =============================================================
%  2) HIGH-PASS FILTER
% ==============================================================
Fp_hpf = 2500;      % Passband edge (Hz)
Fs_hpf = 1800;      % Stopband edge (Hz)

Wp_hpf = Fp_hpf / (Fs/2);
Ws_hpf = Fs_hpf / (Fs/2);

[n_hpf, Wn_hpf] = ellipord(Wp_hpf, Ws_hpf, Rp, Rs);
[b_hpf, a_hpf] = ellip(n_hpf, Rp, Rs, Wn_hpf, 'high');

y_hpf = filter(b_hpf, a_hpf, x);
y_hpf = y_hpf / max(abs(y_hpf) + eps);

% Frequency response
[H_hpf, f_hpf] = freqz(b_hpf, a_hpf, 2048, Fs);

% Spectra
Y_hpfmag = abs(fft(y_hpf));

fig2 = figure;
subplot(3,1,1);
plot(t, x, 'b'); hold on;
plot(t, y_hpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('High-Pass Filter: Original vs Filtered (Time Domain)');
legend('Original','Filtered');
grid on;

subplot(3,1,2);
plot(f_hpf, abs(H_hpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('High-Pass Filter Frequency Response');
grid on;

subplot(3,1,3);
plot(fAxis, X_mag, 'b'); hold on;
plot(fAxis, Y_hpfmag, 'r');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('High-Pass Filter: Original vs Filtered (Frequency Domain)');
legend('Original','Filtered');
grid on;

saveas(fig2, 'assignment5_hpf_analysis.png');

%% =============================================================
%  3) BAND-PASS FILTER
% ==============================================================
Fp1_bpf = 500;      % Lower passband edge (Hz)
Fp2_bpf = 2500;     % Upper passband edge (Hz)
Fs1_bpf = 300;      % Lower stopband edge (Hz)
Fs2_bpf = 3200;     % Upper stopband edge (Hz)

Wp_bpf = [Fp1_bpf Fp2_bpf] / (Fs/2);
Ws_bpf = [Fs1_bpf Fs2_bpf] / (Fs/2);

[n_bpf, Wn_bpf] = ellipord(Wp_bpf, Ws_bpf, Rp, Rs);
[b_bpf, a_bpf] = ellip(n_bpf, Rp, Rs, Wn_bpf, 'bandpass');

y_bpf = filter(b_bpf, a_bpf, x);
y_bpf = y_bpf / max(abs(y_bpf) + eps);

% Frequency response
[H_bpf, f_bpf] = freqz(b_bpf, a_bpf, 2048, Fs);

% Spectra
Y_bpfmag = abs(fft(y_bpf));

fig3 = figure;
subplot(3,1,1);
plot(t, x, 'b'); hold on;
plot(t, y_bpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Band-Pass Filter: Original vs Filtered (Time Domain)');
legend('Original','Filtered');
grid on;

subplot(3,1,2);
plot(f_bpf, abs(H_bpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Band-Pass Filter Frequency Response');
grid on;

subplot(3,1,3);
plot(fAxis, X_mag, 'b'); hold on;
plot(fAxis, Y_bpfmag, 'r');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Band-Pass Filter: Original vs Filtered (Frequency Domain)');
legend('Original','Filtered');
grid on;

saveas(fig3, 'assignment5_bpf_analysis.png');

%% =============================================================
%  4) COMB FILTER
% ==============================================================
delay_ms = 20;                  % Delay in milliseconds
gain_cf  = 0.7;                 % Feedforward comb gain

delay_samples = round(Fs * delay_ms / 1000);

b_cf = [1 zeros(1, delay_samples) gain_cf];
a_cf = 1;

y_cf = filter(b_cf, a_cf, x);
y_cf = y_cf / max(abs(y_cf) + eps);

% Frequency response
[H_cf, f_cf] = freqz(b_cf, a_cf, 2048, Fs);

% Spectra
Y_cfmag = abs(fft(y_cf));

fig4 = figure;
subplot(3,1,1);
plot(t, x, 'b'); hold on;
plot(t, y_cf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Comb Filter: Original vs Filtered (Time Domain)');
legend('Original','Filtered');
grid on;

subplot(3,1,2);
plot(f_cf, abs(H_cf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Comb Filter Frequency Response');
grid on;

subplot(3,1,3);
plot(fAxis, X_mag, 'b'); hold on;
plot(fAxis, Y_cfmag, 'r');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Comb Filter: Original vs Filtered (Frequency Domain)');
legend('Original','Filtered');
grid on;

saveas(fig4, 'assignment5_comb_analysis.png');

%% -------------------------------------------------------------
%  Save filtered audio outputs
% --------------------------------------------------------------
audiowrite('assignment5_original_segment.wav', x, Fs);
audiowrite('assignment5_lpf_output.wav', y_lpf, Fs);
audiowrite('assignment5_hpf_output.wav', y_hpf, Fs);
audiowrite('assignment5_bpf_output.wav', y_bpf, Fs);
audiowrite('assignment5_comb_output.wav', y_cf, Fs);

%% -------------------------------------------------------------
%  Playback sequence
% --------------------------------------------------------------
disp('Playing original segment...');
soundsc(x, Fs);
pause(length(x)/Fs + 1);

disp('Playing LPF output...');
soundsc(y_lpf, Fs);
pause(length(y_lpf)/Fs + 1);

disp('Playing HPF output...');
soundsc(y_hpf, Fs);
pause(length(y_hpf)/Fs + 1);

disp('Playing BPF output...');
soundsc(y_bpf, Fs);
pause(length(y_bpf)/Fs + 1);

disp('Playing Comb filter output...');
soundsc(y_cf, Fs);

disp('Assignment 5 processing complete.');
