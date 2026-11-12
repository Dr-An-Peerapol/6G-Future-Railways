
%% SIMO Lesson 1: Understanding the Steering Vector

clear; clc;

N = 4;                % Number of receive antennas
theta = 30;           % Arrival angle in degrees
d_lambda = 0.5;       % Antenna spacing / wavelength

n = (0:N-1).';

% Spatial phase difference
psi = 2*pi*d_lambda*sind(theta);

% Array steering vector
a = exp(-1j*n*psi);

disp('Steering vector:');
disp(a);

disp('Array power gain a^H*a:');
disp(real(a'*a));

% Visualize phases across antennas
figure('Color','w');

stem(1:N, angle(a)*180/pi, 'filled', ...
    'LineWidth',1.5);
grid on;
xlabel('Antenna number');
ylabel('Received phase (degrees)');
title('Spatial Phase Signature of an Incoming Wave');
xticks(1:N);
