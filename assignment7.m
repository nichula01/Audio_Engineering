% ===============================================================
% EE599 - Audio Engineering and Acoustics
% Assignment 07
% Form Filters Using Z-plane Poles and Zeros
%
% Student Name : Nichula Wasalathilaka
% Reg. No.     : E/20/425
%
% Instrument   : Guitar
%
% Description:
% This script forms four filters directly using z-plane pole-zero
% placement:
%   1) Low-Pass Filter (LPF)
%   2) High-Pass Filter (HPF)
%   3) Band-Pass Filter (BPF)
%   4) Comb Filter
%
% For each filter, it generates one figure containing:
%   - time-domain comparison
%   - magnitude response
%   - frequency-domain comparison
%   - pole-zero plot
%
% All figures are saved as PNG files.
% ===============================================================

clc;
clear;
close all;

% Safer rendering for MATLAB Online
set(groot, 'defaultFigureRenderer', 'painters');
set(groot, 'defaultFigureVisible', 'off');

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
[fileData, fs] = audioread(audioFileName);

if size(fileData,2) > 1
    x = mean(fileData,2);
else
    x = fileData;
end

x = x(:);
x = x / max(abs(x) + eps);

start_sec = 0;
dur_sec   = 4;   % keep shorter for faster saving

start_sample = floor(start_sec * fs) + 1;
end_sample   = min(length(x), floor((start_sec + dur_sec) * fs));

x = x(start_sample:end_sample);
x = x / max(abs(x) + eps);

t = (0:length(x)-1)/fs;

X = fft(x);
f_axis = (0:length(X)-1)*(fs/length(X));

%% =============================================================
%  1) LOW-PASS FILTER USING Z-PLANE
% =============================================================
p_lpf = [ ...
    0.90*exp(1j*0.10*pi), 0.90*exp(-1j*0.10*pi), ...
    0.85*exp(1j*0.20*pi), 0.85*exp(-1j*0.20*pi) ...
].';

z_lpf = [-1; -1; -1; -1];

[b_lpf, a_lpf] = zp2tf(z_lpf, p_lpf, 1);
y_lpf = filter(b_lpf, a_lpf, x);
y_lpf = y_lpf / max(abs(y_lpf) + eps);

[h_lpf, w_lpf] = freqz(b_lpf, a_lpf, 1024);
f_lpf = w_lpf * fs / (2*pi);

Y_lpf = fft(y_lpf);

fig1 = figure('Visible','off','Renderer','painters','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t, x, 'b'); hold on;
plot(t, y_lpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('LPF: Original vs Filtered');
legend('Original','Filtered');
grid on;
xlim([0 min(0.1, t(end))]);

subplot(2,2,2);
plot(f_lpf, abs(h_lpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('LPF Magnitude Response');
grid on;

subplot(2,2,3);
plot(f_axis, abs(X), 'b'); hold on;
plot(f_axis, abs(Y_lpf), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('LPF Frequency Domain');
legend('Original','Filtered');
xlim([0 fs/2]);
grid on;

subplot(2,2,4);
zplane(b_lpf, a_lpf);
title('LPF Pole-Zero Plot');

print(fig1, 'assignment7_lpf_zplane.png', '-dpng', '-r200');
close(fig1);

%% =============================================================
%  2) HIGH-PASS FILTER USING Z-PLANE
% =============================================================
p_hpf = [ ...
    0.90*exp(1j*0.90*pi), 0.90*exp(-1j*0.90*pi), ...
    0.85*exp(1j*0.80*pi), 0.85*exp(-1j*0.80*pi) ...
].';

z_hpf = [1; 1; 1; 1];

[b_hpf, a_hpf] = zp2tf(z_hpf, p_hpf, 1);
y_hpf = filter(b_hpf, a_hpf, x);
y_hpf = y_hpf / max(abs(y_hpf) + eps);

[h_hpf, w_hpf] = freqz(b_hpf, a_hpf, 1024);
f_hpf = w_hpf * fs / (2*pi);

Y_hpf = fft(y_hpf);

fig2 = figure('Visible','off','Renderer','painters','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t, x, 'b'); hold on;
plot(t, y_hpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('HPF: Original vs Filtered');
legend('Original','Filtered');
grid on;
xlim([0 min(0.1, t(end))]);

subplot(2,2,2);
plot(f_hpf, abs(h_hpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('HPF Magnitude Response');
grid on;

subplot(2,2,3);
plot(f_axis, abs(X), 'b'); hold on;
plot(f_axis, abs(Y_hpf), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('HPF Frequency Domain');
legend('Original','Filtered');
xlim([0 fs/2]);
grid on;

subplot(2,2,4);
zplane(b_hpf, a_hpf);
title('HPF Pole-Zero Plot');

print(fig2, 'assignment7_hpf_zplane.png', '-dpng', '-r200');
close(fig2);

%% =============================================================
%  3) BAND-PASS FILTER USING Z-PLANE
% =============================================================
theta = 0.40*pi;

p_bpf = [ ...
    0.90*exp(1j*theta),        0.90*exp(-1j*theta), ...
    0.85*exp(1j*(theta+0.10)), 0.85*exp(-1j*(theta+0.10)) ...
].';

z_bpf = [1; -1; 1; -1];

[b_bpf, a_bpf] = zp2tf(z_bpf, p_bpf, 1);
y_bpf = filter(b_bpf, a_bpf, x);
y_bpf = y_bpf / max(abs(y_bpf) + eps);

[h_bpf, w_bpf] = freqz(b_bpf, a_bpf, 1024);
f_bpf = w_bpf * fs / (2*pi);

Y_bpf = fft(y_bpf);

fig3 = figure('Visible','off','Renderer','painters','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t, x, 'b'); hold on;
plot(t, y_bpf, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('BPF: Original vs Filtered');
legend('Original','Filtered');
grid on;
xlim([0 min(0.1, t(end))]);

subplot(2,2,2);
plot(f_bpf, abs(h_bpf), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('BPF Magnitude Response');
grid on;

subplot(2,2,3);
plot(f_axis, abs(X), 'b'); hold on;
plot(f_axis, abs(Y_bpf), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('BPF Frequency Domain');
legend('Original','Filtered');
xlim([0 fs/2]);
grid on;

subplot(2,2,4);
zplane(b_bpf, a_bpf);
title('BPF Pole-Zero Plot');

print(fig3, 'assignment7_bpf_zplane.png', '-dpng', '-r200');
close(fig3);

%% =============================================================
%  4) COMB FILTER USING Z-PLANE
% =============================================================
delay_ms = 20;
N_delay  = round(0.02 * fs);
alpha_ff = 0.8;
alpha_fb = 0.7;

% Feedforward comb for clearer pole-zero interpretation
b_comb = [1 zeros(1, N_delay-1) -alpha_ff];
a_comb = 1;

y_comb = filter(b_comb, a_comb, x);
y_comb = y_comb / max(abs(y_comb) + eps);

[h_comb, w_comb] = freqz(b_comb, a_comb, 1024);
f_comb = w_comb * fs / (2*pi);

Y_comb = fft(y_comb);

fig4 = figure('Visible','off','Renderer','painters','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(t, x, 'b'); hold on;
plot(t, y_comb, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('Comb Filter: Original vs Filtered');
legend('Original','Filtered');
grid on;
xlim([0 min(0.1, t(end))]);

subplot(2,2,2);
plot(f_comb, abs(h_comb), 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Comb Filter Magnitude Response');
grid on;

subplot(2,2,3);
plot(f_axis, abs(X), 'b'); hold on;
plot(f_axis, abs(Y_comb), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Comb Filter Frequency Domain');
legend('Original','Filtered');
xlim([0 fs/2]);
grid on;

subplot(2,2,4);
zplane(b_comb, a_comb);
title('Comb Filter Pole-Zero Plot');

print(fig4, 'assignment7_comb_zplane.png', '-dpng', '-r200');
close(fig4);

%% -------------------------------------------------------------
%  Save audio outputs
% --------------------------------------------------------------
audiowrite('assignment7_original_segment.wav', x, fs);
audiowrite('assignment7_lpf_output.wav', y_lpf, fs);
audiowrite('assignment7_hpf_output.wav', y_hpf, fs);
audiowrite('assignment7_bpf_output.wav', y_bpf, fs);
audiowrite('assignment7_comb_output.wav', y_comb, fs);

disp('Assignment 7 processing complete.');
disp('Saved figures:');
disp('  assignment7_lpf_zplane.png');
disp('  assignment7_hpf_zplane.png');
disp('  assignment7_bpf_zplane.png');
disp('  assignment7_comb_zplane.png');
