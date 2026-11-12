
%% SIMO Lesson 5:
% Four-element linear vs 2x2 planar array

clear; clc;

d = 0.5;  % Spacing normalized by wavelength

% 1x4 linear array positions (x,y)
posULA = d*[0 0; 1 0; 2 0; 3 0];

% 2x2 planar array positions (x,y)
posUPA = d*[0 0; 1 0; 0 1; 1 1];

% General steering vector for arbitrary array positions
steer = @(pos,theta,phi) exp( ...
    -1j*2*pi*( ...
    pos(:,1)*sind(theta)*cosd(phi) + ...
    pos(:,2)*sind(theta)*sind(phi)));

rho = @(a,b) abs(a'*b)/(norm(a)*norm(b));

% Two example geometries:
% [theta1 phi1 theta2 phi2]
cases = [45 45 45 -45;
         30  0 30  90];

for k = 1:size(cases,1)

    s = cases(k,:);

    % Linear array
    a1L = steer(posULA,s(1),s(2));
    a2L = steer(posULA,s(3),s(4));

    % Planar array
    a1P = steer(posUPA,s(1),s(2));
    a2P = steer(posUPA,s(3),s(4));

    AL = [a1L a2L];
    AP = [a1P a2P];

    GL = (AL'*AL)/4;
    GP = (AP'*AP)/4;

    fprintf('\nCase %d\n',k);

    fprintf('Linear correlation: %.4f\n', ...
        rho(a1L,a2L));

    fprintf('Planar correlation: %.4f\n', ...
        rho(a1P,a2P));

    fprintf('Linear Gram eigenvalues:\n');
    disp(sort(real(eig(GL))).');

    fprintf('Planar Gram eigenvalues:\n');
    disp(sort(real(eig(GP))).');

end
