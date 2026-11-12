
%% SIMO Lesson 6: Spatial Channel Covariance

clear; clc;

N = 4;
d = 0.5;

posULA = d*[0 0; 1 0; 2 0; 3 0];
posUPA = d*[0 0; 1 0; 0 1; 1 1];

steer = @(pos,theta,phi) exp(-1j*2*pi*( ...
    pos(:,1)*sind(theta)*cosd(phi) + ...
    pos(:,2)*sind(theta)*sind(phi)));

% Same propagation directions as Case 1
theta1 = 45; phi1 = 45;
theta2 = 45; phi2 = -45;

% Two equal-power uncorrelated paths
Rs = 0.5*eye(2);

% Steering matrices
Aula = [steer(posULA,theta1,phi1), ...
        steer(posULA,theta2,phi2)];

Aupa = [steer(posUPA,theta1,phi1), ...
        steer(posUPA,theta2,phi2)];

% Spatial channel covariance matrices
Rula = Aula*Rs*Aula';
Rupa = Aupa*Rs*Aupa';

% Eigenvalue decompositions
[Uula,Dula] = eig(Rula);
[Uupa,Dupa] = eig(Rupa);

fprintf('Linear covariance eigenvalues:\n');
disp(sort(real(diag(Dula)),'descend').');

fprintf('Planar covariance eigenvalues:\n');
disp(sort(real(diag(Dupa)),'descend').');

fprintf('Total power, linear = %.4f\n',real(trace(Rula)));
fprintf('Total power, planar = %.4f\n',real(trace(Rupa)));

% Compare eigenvalue spectra
figure('Color','w');
plot(1:4,sort(real(eig(Rula)),'descend'), ...
    'o-','LineWidth',2);
hold on;
plot(1:4,sort(real(eig(Rupa)),'descend'), ...
    's-','LineWidth',2);
grid on;
xlabel('Eigenvalue index');
ylabel('Spatial covariance eigenvalue');
title('Spatial Eigenmodes: Linear vs Planar Array');
legend('1x4 Linear','2x2 Planar','Location','best');
xticks(1:4);
