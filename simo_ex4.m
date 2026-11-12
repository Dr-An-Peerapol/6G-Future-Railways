
%% SIMO Lesson 4: Gram Matrix of Steering Vectors

clear; clc;

N = 4;
n = (0:N-1).';

a = @(theta) exp(-1j*pi*n*sind(theta));

% Case 1: Orthogonal directions
A1 = [a(30), a(-30)];
G1 = (A1'*A1)/N;

% Case 2: Closely spaced directions
A2 = [a(30), a(35)];
G2 = (A2'*A2)/N;

disp('Case 1: |G| for 30 and -30 degrees');
disp(abs(G1));

disp('Case 2: |G| for 30 and 35 degrees');
disp(abs(G2));

disp('Eigenvalues of Gram matrix, Case 1');
disp(eig(G1));

disp('Eigenvalues of Gram matrix, Case 2');
disp(eig(G2));
