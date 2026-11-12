
%% SIMO Lesson 3: Steering Vector Correlation

clear; close all; clc;

N = 4;
d_lambda = 0.5;

n = (0:N-1).';

% First wave arrives from 30 degrees
theta1 = 30;

a1 = exp(-1j*2*pi*d_lambda*n*sind(theta1));

% Test other arrival angles
theta2 = -90:0.5:90;

rho = zeros(size(theta2));

for k = 1:length(theta2)

    a2 = exp(-1j*2*pi*d_lambda*n*sind(theta2(k)));

    % Normalized inner product
    rho(k) = abs(a1'*a2)/(norm(a1)*norm(a2));

end

figure('Color','w');

plot(theta2,rho,'b-','LineWidth',2);
grid on;

xlabel('Second wave arrival angle (degrees)');
ylabel('Normalized spatial correlation');

title('Spatial Correlation Between Two Steering Vectors');

xline(theta1,'r--','Reference angle');

ylim([0 1.05]);

% Examine three specific angles
testAngles = [-30 30 35];

fprintf('\nSpatial correlation with theta1 = 30 deg\n');

for theta = testAngles

    a2 = exp(-1j*2*pi*d_lambda*n*sind(theta));

    correlation = abs(a1'*a2)/(norm(a1)*norm(a2));

    fprintf('theta2 = %+5.1f deg : rho = %.4f\n', ...
        theta, correlation);

end
