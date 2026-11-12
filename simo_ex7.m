
%% EXPERIMENT 7: SIMO 1x4 WITH TWO COMMUNICATION LINKS
% Link A: Two propagation rays, 1000 Hz Doppler separation
% Link B: One constant-power propagation ray
%
% Compare:
%   - Single antenna
%   - Perfect MRC
%   - Stale MRC
%   - Predictive MRC with Doppler errors

clear; close all; clc;

%% 1. Railway parameters

c = 299792458;
fc = 28e9;
v = 350/3.6;

lambda = c/fc;
fDmax = v/lambda;

fprintf('Maximum Doppler = %.1f Hz\n',fDmax);

%% 2. Four-element linear array

N = 4;
dLambda = 0.5;

n = (0:N-1).';

steer = @(theta) exp( ...
    -1j*2*pi*dLambda*n*sind(theta));

% Arrival directions
a1 = steer(+30);
a2 = steer(-30);
aB = steer(+10);

fprintf('Spatial correlation of A rays = %.5f\n', ...
    abs(a1'*a2)/(norm(a1)*norm(a2)));

%% 3. Doppler frequencies

f1 = 0.9*fDmax;
f2 = f1-1000;
fB = f1-500;

% Imperfect Doppler estimates used for prediction
f1Pred = f1+50;
f2Pred = f2-50;

%% 4. Prediction horizons

timeMs = 0:0.01:3;
t = timeMs*1e-3;

%% 5. Actual channels: 4 antennas x time

hA = ( ...
    a1*exp(1j*2*pi*f1*t) + ...
    a2*exp(1j*2*pi*f2*t))/sqrt(2);

hB = sqrt(0.8)*aB*exp(1j*2*pi*fB*t);

%% 6. Predicted Link A channel

hAPred = ( ...
    a1*exp(1j*2*pi*f1Pred*t) + ...
    a2*exp(1j*2*pi*f2Pred*t))/sqrt(2);

%% 7. Single-antenna power

gainA1 = abs(hA(1,:)).^2;
gainB1 = abs(hB(1,:)).^2;

%% 8. Perfect four-antenna MRC

gainAMRC = sum(abs(hA).^2,1);
gainBMRC = sum(abs(hB).^2,1);

%% 9. Stale MRC: use channel at time zero

wStale = hA(:,1);

gainAStale = abs(wStale'*hA).^2 / ...
    (wStale'*wStale);

gainAStale = real(gainAStale);

%% 10. Predictive MRC

% At each horizon, use its predicted channel vector.
% The prediction is imperfect due to Doppler errors.

numerator = abs(sum(conj(hAPred).*hA,1)).^2;

denominator = sum(abs(hAPred).^2,1);

gainAPred = numerator./denominator;

%% 11. Compare link selections

% Perfect information about single-antenna powers
chooseSISO = ones(size(t));
chooseSISO(gainB1>gainA1) = 2;

% Perfect information about four-antenna MRC powers
chooseSIMO = ones(size(t));
chooseSIMO(gainBMRC>gainAMRC) = 2;

%% 12. Display figures

figure('Color','w','Position',[100 100 1150 750]);

subplot(2,2,1);

plot(timeMs,abs(hA).^2,'LineWidth',1.5);

grid on;
xlabel('Time (ms)');
ylabel('Received power');
title('Link A: Power at Each Antenna');
legend('Rx1','Rx2','Rx3','Rx4', ...
    'Location','best');

subplot(2,2,2);

plot(timeMs,gainA1,'b-','LineWidth',1.7);
hold on;
plot(timeMs,gainB1,'r--','LineWidth',1.7);

plot(timeMs,gainAMRC,'k-','LineWidth',2);
plot(timeMs,gainBMRC,'m--','LineWidth',2);

grid on;
xlabel('Time (ms)');
ylabel('Effective channel power gain');
title('One Antenna vs Four-Antenna MRC');

legend('Link A: 1 antenna', ...
       'Link B: 1 antenna', ...
       'Link A: 4 antennas', ...
       'Link B: 4 antennas', ...
       'Location','best');

subplot(2,2,3);

plot(timeMs,gainAMRC,'k-','LineWidth',2);
hold on;
plot(timeMs,gainAStale,'r--','LineWidth',1.7);
plot(timeMs,gainAPred,'b-','LineWidth',1.7);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Effective channel gain');
title('Link A: Effect of Aging CSI');

legend('Perfect MRC', ...
       'Stale MRC', ...
       'Predictive MRC', ...
       'Location','best');

subplot(2,2,4);

stairs(timeMs,chooseSISO,'r-','LineWidth',1.5);
hold on;
stairs(timeMs,chooseSIMO,'b--','LineWidth',1.5);

grid on;
ylim([0.7 2.3]);

yticks([1 2]);
yticklabels({'Link A','Link B'});

xlabel('Time (ms)');
title('Ideal Link Selection: SISO vs SIMO');

legend('Single antenna', ...
       'Four-antenna MRC', ...
       'Location','best');

sgtitle('SIMO Reception with Two-Ray Multipath');

%% 13. Examine the deep fade at 0.5 ms

[~,idx] = min(abs(timeMs-0.5));

fprintf('\n=== AT 0.5 ms ===\n');

fprintf('Link A, antenna 1 power : %.5f\n', ...
    gainA1(idx));

fprintf('Link A, perfect MRC gain: %.5f\n', ...
    gainAMRC(idx));

fprintf('Link A, stale MRC gain  : %.5f\n', ...
    gainAStale(idx));

fprintf('Link A, predicted gain  : %.5f\n', ...
    gainAPred(idx));

fprintf('Link B, perfect MRC gain: %.5f\n', ...
    gainBMRC(idx));

%% 14. QPSK reception test at 0.5 ms

% Use 0 dB per-antenna reference SNR so that
% the improvement is visible with finite symbols.

snrQPSK = 1;     % 0 dB
Nsym = 50000;

bitsI = randi([0 1],1,Nsym);
bitsQ = randi([0 1],1,Nsym);

x = ((2*bitsI-1)+1j*(2*bitsQ-1))/sqrt(2);

% Independent complex Gaussian noise on four branches
noise = sqrt(1/(2*snrQPSK))* ...
    (randn(N,Nsym)+1j*randn(N,Nsym));

% Receive through Link A at t = 0.5 ms
hNow = hA(:,idx);

y = hNow*x+noise;

% Single-antenna coherent detection
xHat1 = y(1,:)/hNow(1);

% Four-antenna perfect MRC
xHat4 = (hNow'*y)/(hNow'*hNow);

% QPSK sign decisions
bI1 = real(xHat1)>0;
bQ1 = imag(xHat1)>0;

bI4 = real(xHat4)>0;
bQ4 = imag(xHat4)>0;

BER1 = (sum(bI1~=bitsI)+sum(bQ1~=bitsQ))/(2*Nsym);
BER4 = (sum(bI4~=bitsI)+sum(bQ4~=bitsQ))/(2*Nsym);

fprintf('\n=== QPSK AT 0.5 ms, 0 dB SNR ===\n');
fprintf('Single-antenna BER: %.5f\n',BER1);
fprintf('Four-antenna BER  : %.5f\n',BER4);

fprintf('\nExperiment complete.\n');
