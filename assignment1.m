% ---------------------------------------------------------
% EE599 Assignment 3
% Instrument 1 - Harmonic string-like resonator
% Created by: [Your Name]
%
% This script creates a simple instrument-like sound by
% placing complex poles at a fundamental frequency and its
% harmonics, then shaping the tone using zeros.
% ---------------------------------------------------------

clc;
clear;
close all;

%% Basic settings
Fs = 8000;              % Sampling frequency in Hz
N  = 20000;             % Length of impulse response
f0 = 330;               % Fundamental frequency in Hz

%% Pole radii
% Radii are chosen close to 1 to create a sustained ringing sound.
r1 = 0.996;
r2 = 0.994;
r3 = 0.992;

%% Convert frequencies to digital angular frequencies
w1 = 2*pi*f0/Fs;        % Fundamental
w2 = 2*pi*(2*f0)/Fs;    % 2nd harmonic
w3 = 2*pi*(3*f0)/Fs;    % 3rd harmonic

%% Define poles
% Complex-conjugate pole pairs are used for real-valued output.
poles = [ ...
    r1*exp(1j*w1), r1*exp(-1j*w1), ...
    r2*exp(1j*w2), r2*exp(-1j*w2), ...
    r3*exp(1j*w3), r3*exp(-1j*w3) ...
];

%% Define zeros
% Zeros are used to shape the tone and reduce harsh high-frequency content.
zeros_z = [-0.75, -0.75];

%% Convert pole-zero locations to filter coefficients
a = poly(poles);        % Denominator coefficients
b = poly(zeros_z);      % Numerator coefficients

%% Generate a short excitation signal
% A short burst of random noise makes the result more plucked and natural.
x = 0.02 * randn(N,1);
x(250:end) = 0;

%% Filter the excitation through the resonator
y = filter(b, a, x);

%% Normalize output to avoid clipping
y = y / max(abs(y));

%% Time-domain plot
figure;
plot(y, 'LineWidth', 1);
xlabel('Sample Index');
ylabel('Amplitude');
title('Instrument 1 - Time Domain Response');
grid on;

%% Pole-zero plot
figure;
zplane(b, a);
title('Instrument 1 - Pole-Zero Plot');

%% Play the sound
soundsc(y, Fs);
