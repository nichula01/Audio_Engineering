% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 06
% Form an Equalizer Combining LPF, HPF and BPF
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script builds a simple 3-band equalizer using:
%   1) Low-Pass Filter  -> bass band
%   2) Band-Pass Filter -> mid band
%   3) High-Pass Filter -> treble band
%
% The filtered bands are scaled using separate gains and then
% combined to form the equalized output signal.
%
% Outputs:
%   - Time-domain comparison figure
%   - Equalizer frequency response figure
%   - Frequency-domain comparison figure
%   - Original and equalized WAV files
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

N = length(x);
t = (0:N-1) / Fs;

%% -------------------------------------------------------------
%  Filter design settings
% --------------------------------------------------------------
Rp = 1;      % Passband ripple in dB
Rs = 40;     % Stopband attenuation in dB

%% -------------------------------------------------------------
%  Low-Pass Filter (Bass band)
% --------------------------------------------------------------
Wp_lp = 800 / (Fs/2);
Ws_lp = 1000 / (Fs/2);

[n_lp, Wn_lp] = ellipord(Wp_lp, Ws_lp, Rp, Rs);
[b_lp, a_lp] = ellip(n_lp, Rp, Rs, Wn_lp, 'low');

%% -------------------------------------------------------------
%  Band-Pass Filter (Mid band)
% --------------------------------------------------------------
Wp_bp = [1000 3000] / (Fs/2);
Ws_bp = [800 3500] / (Fs/2);

[n_bp, Wn_bp] = ellipord(Wp_bp, Ws_bp, Rp, Rs);
[b_bp, a_bp] = ellip(n_bp, Rp, Rs, Wn_bp, 'bandpass');

%% -------------------------------------------------------------
%  High-Pass Filter (Treble band)
% --------------------------------------------------------------
Wp_hp = 3000 / (Fs/2);
Ws_hp = 2500 / (Fs/2);

[n_hp, Wn_hp] = ellipord(Wp_hp, Ws_hp, Rp, Rs);
[b_hp, a_hp] = ellip(n_hp, Rp, Rs, Wn_hp, 'high');

%% -------------------------------------------------------------
%  Filter the signal into three bands
% --------------------------------------------------------------
y_lp = filter(b_lp, a_lp, x);
y_bp = filter(b_bp, a_bp, x);
y_hp = filter(b_hp, a_hp, x);

%% -------------------------------------------------------------
%  Equalizer gains
% --------------------------------------------------------------
G_lp = 1.2;     % Boost bass
G_bp = 1.0;     % Keep mid band unchanged
G_hp = 0.8;     % Slightly reduce treble

%% -------------------------------------------------------------
%  Combine the bands to form equalized output
% --------------------------------------------------------------
y_eq = G_lp * y_lp + G_bp * y_bp + G_hp * y_hp;
y_eq = y_eq / max(abs(y_eq) + eps);

%% -------------------------------------------------------------
%  Plot 1: Time-domain comparison
% --------------------------------------------------------------
fig1 = figure;
plot(t, x, 'b'); hold on;
plot(t, y_eq, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Original vs Equalized Signal (Time Domain)');
legend('Original', 'Equalized');
grid on;
xlim([0 min(0.1, t(end))]);   % Zoom slightly for clarity

saveas(fig1, 'assignment6_time_comparison.png');

%% -------------------------------------------------------------
%  Plot 2: Equalizer frequency response
% --------------------------------------------------------------
[H_lp, f] = freqz(b_lp, a_lp, 2048, Fs);
H_bp      = freqz(b_bp, a_bp, 2048, Fs);
H_hp      = freqz(b_hp, a_hp, 2048, Fs);

H_eq = G_lp * H_lp + G_bp * H_bp + G_hp * H_hp;

fig2 = figure;
plot(f, 20*log10(abs(H_lp) + 1e-6), 'LineWidth', 1.1); hold on;
plot(f, 20*log10(abs(H_bp) + 1e-6), 'LineWidth', 1.1);
plot(f, 20*log10(abs(H_hp) + 1e-6), 'LineWidth', 1.1);
plot(f, 20*log10(abs(H_eq) + 1e-6), 'k', 'LineWidth', 1.6);

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Equalizer Frequency Response');
legend('LPF (Bass)', 'BPF (Mid)', 'HPF (Treble)', 'Combined EQ');
grid on;

saveas(fig2, 'assignment6_eq_response.png');

%% -------------------------------------------------------------
%  Plot 3: Frequency-domain comparison
% --------------------------------------------------------------
X = abs(fft(x));
Y = abs(fft(y_eq));
f_fft = (0:N-1) * (Fs/N);

fig3 = figure;
plot(f_fft, X, 'b'); hold on;
plot(f_fft, Y, 'r');
xlim([0 Fs/2]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Original vs Equalized Signal (Frequency Domain)');
legend('Original', 'Equalized');
grid on;

saveas(fig3, 'assignment6_frequency_comparison.png');

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment6_original_segment.wav', x, Fs);
audiowrite('assignment6_equalized_output.wav', y_eq, Fs);

%% -------------------------------------------------------------
%  Playback
% --------------------------------------------------------------
disp('Playing original segment...');
soundsc(x, Fs);
pause(length(x)/Fs + 1);

disp('Playing equalized output...');
soundsc(y_eq, Fs);

disp('Assignment 6 processing complete.');
