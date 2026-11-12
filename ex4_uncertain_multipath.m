
%% EXPERIMENT 4: TWO UNCERTAIN MULTIPATH LINKS
% High-speed railway, 350 km/h, 28 GHz
%
% Link A: 5 propagation paths
% Link B: 4 propagation paths
%
% Both links have:
%   - Different path powers and Doppler shifts
%   - Uncertain Doppler
%   - Multipath fading
%   - Possible LoS blockage
%
% Compare:
%   1. Stale link selection
%   2. Predictive rate-based selection
%   3. Risk-aware selection
%   4. Fresh-pilot-assisted selection
%   5. Perfect future CSI (oracle)
%
% Base MATLAB only. MATLAB R2016b or newer.

clear; close all; clc;
rng(42);

%% 1. Physical parameters

c = 299792458;
fc = 28e9;
v = 350/3.6;
lambda = c/fc;
fDmax = v/lambda;

fprintf('Train speed: %.2f km/h\n',v*3.6);
fprintf('Maximum Doppler: %.2f Hz\n',fDmax);

%% 2. Two different multipath links

% Arrival angles relative to train velocity (degrees)
thetaA = [10 48 95 135 172];
thetaB = [25 70 115 155];

% Path powers: LoS followed by reflected paths
powerA = [0.50 0.18 0.15 0.10 0.07];
powerB = [0.45 0.27 0.18 0.10];

% Initial phases in radians
phaseA = [0 0.35 -0.6 1.2 -1.7];
phaseB = [0.2 1.1 -0.9 2.0];

% Nominal Doppler for each path
fNomA = fDmax*cosd(thetaA);
fNomB = fDmax*cosd(thetaB);

% Complex path gains at time zero
aA = sqrt(powerA(:)).*exp(1j*phaseA(:));
aB = sqrt(powerB(:)).*exp(1j*phaseB(:));

LA = numel(aA);
LB = numel(aB);

hA0 = sum(aA);
hB0 = sum(aB);

fprintf('Initial Link A power: %.3f\n',abs(hA0)^2);
fprintf('Initial Link B power: %.3f\n',abs(hB0)^2);

%% 3. Uncertainty and communication requirements

% Shared train-speed estimation error:
% introduces correlation across paths and links
sigmaSpeedRelative = 0.004;

% Additional path-dependent Doppler uncertainty
sigmaResidualFD = 40;    % Hz

% Reference signal-to-noise ratio
snrLinear = 10;          % 10 dB

% Communication requirement
Rmin = 2.5;              % bit/s/Hz
epsilon = 0.05;          % target outage probability

% Expected-regret tolerance
regretBudget = 0.10;     % bit/s/Hz

% Future pilot measurement settings
pilotNoiseVar = 0.025;
pilotOverhead = 0.05;

%% 4. Monte Carlo experiment

Nmodel = 1500;           % For controller prediction
Ntest = 3000;            % Independent evaluation set
N = Nmodel + Ntest;

modelIdx = 1:Nmodel;
testIdx = Nmodel+(1:Ntest);

% Each scenario keeps the same uncertain Doppler
% parameters throughout the simulated trajectory.

speedError = sigmaSpeedRelative*randn(N,1);

errA = speedError.*fNomA ...
    + sigmaResidualFD*randn(N,LA);

errB = speedError.*fNomB ...
    + sigmaResidualFD*randn(N,LB);

fA = fNomA + errA;
fB = fNomB + errB;

%% 5. Random blockage events

% These are illustrative stochastic blockage events.
% The predictor knows their probability model, not
% which future realization will occur.

eventA = rand(N,1)<0.30;
eventB = rand(N,1)<0.30;
eventCommon = rand(N,1)<0.12;

startA = 1.5e-3 + 3e-3*rand(N,1);
startB = 1.5e-3 + 3e-3*rand(N,1);
startCommon = 2e-3 + 3e-3*rand(N,1);

endA = startA + 1.5e-3 + 1.5e-3*rand(N,1);
endB = startB + 1.5e-3 + 1.5e-3*rand(N,1);
endCommon = startCommon + 2e-3;

% Blockage attenuates the dominant path.
% Scattered paths are still present.

%% 6. Future prediction horizons

timeMs = (0:0.02:10).';
timeSec = timeMs*1e-3;
Nt = numel(timeMs);

% Results
meanPowerA = zeros(Nt,1);
meanPowerB = zeros(Nt,1);
varChanA = zeros(Nt,1);
varChanB = zeros(Nt,1);

meanRateA = zeros(Nt,1);
meanRateB = zeros(Nt,1);
testRateA = zeros(Nt,1);
testRateB = zeros(Nt,1);

pOutA = zeros(Nt,1);
pOutB = zeros(Nt,1);

choicePred = zeros(Nt,1);
choiceRisk = zeros(Nt,1);
bothInfeasible = false(Nt,1);

rateStale = zeros(Nt,1);
ratePred = zeros(Nt,1);
rateRisk = zeros(Nt,1);
ratePilot = zeros(Nt,1);
rateOracle = zeros(Nt,1);

outStale = zeros(Nt,1);
outPred = zeros(Nt,1);
outRisk = zeros(Nt,1);
outPilot = zeros(Nt,1);

% Sample channel power trajectories
Nshow = 8;
sampleA = zeros(Nshow,Nt);
sampleB = zeros(Nshow,Nt);

%% 7. Determine the stale decision

if abs(hA0)^2 >= abs(hB0)^2
    staleLink = 1;
else
    staleLink = 2;
end

fprintf('Initially selected link: %d\n',staleLink);

%% 8. Main simulation

for k = 1:Nt

    t = timeSec(k);

    %% A. Calculate both physical multipath channels

    % Sum of complex Doppler-shifted paths
    expA = exp(1j*2*pi*t*fA);
    expB = exp(1j*2*pi*t*fB);

    hA = expA*aA;
    hB = expB*aB;

    % Extract direct-path contribution
    losA = expA(:,1)*aA(1);
    losB = expB(:,1)*aB(1);

    %% B. Apply time-varying blockage

    blockedA = eventA & (t>=startA) & (t<endA);
    blockedB = eventB & (t>=startB) & (t<endB);

    blockedCommon = eventCommon & ...
        (t>=startCommon) & (t<endCommon);

    % Amplitude attenuation of LoS contribution
    attA = (1-0.85*double(blockedA)) .* ...
           (1-0.40*double(blockedCommon));

    attB = (1-0.80*double(blockedB)) .* ...
           (1-0.55*double(blockedCommon));

    hA = hA+(attA-1).*losA;
    hB = hB+(attB-1).*losB;

    %% C. Separate prediction and evaluation samples

    hAm = hA(modelIdx);
    hBm = hB(modelIdx);

    hAt = hA(testIdx);
    hBt = hB(testIdx);

    PAmodel = abs(hAm).^2;
    PBmodel = abs(hBm).^2;

    PAtest = abs(hAt).^2;
    PBtest = abs(hBt).^2;

    sampleA(:,k) = PAtest(1:Nshow);
    sampleB(:,k) = PBtest(1:Nshow);

    %% D. Future channel moments

    meanPowerA(k) = mean(PAmodel);
    meanPowerB(k) = mean(PBmodel);

    muA = mean(hAm);
    muB = mean(hBm);

    varChanA(k) = mean(abs(hAm-muA).^2);
    varChanB(k) = mean(abs(hBm-muB).^2);

    %% E. Predict spectral efficiency and outage

    RAm = log2(1+snrLinear*PAmodel);
    RBm = log2(1+snrLinear*PBmodel);

    RAt = log2(1+snrLinear*PAtest);
    RBt = log2(1+snrLinear*PBtest);

    meanRateA(k) = mean(RAm);
    meanRateB(k) = mean(RBm);

    testRateA(k) = mean(RAt);
    testRateB(k) = mean(RBt);

    pOutA(k) = mean(RAm<Rmin);
    pOutB(k) = mean(RBm<Rmin);

    %% F. Stale controller

    if staleLink==1
        rateStale(k) = mean(RAt);
        outStale(k) = mean(RAt<Rmin);
    else
        rateStale(k) = mean(RBt);
        outStale(k) = mean(RBt<Rmin);
    end

    %% G. Predictive controller: highest expected rate

    if meanRateA(k)>=meanRateB(k)
        choicePred(k) = 1;
        ratePred(k) = mean(RAt);
        outPred(k) = mean(RAt<Rmin);
    else
        choicePred(k) = 2;
        ratePred(k) = mean(RBt);
        outPred(k) = mean(RBt<Rmin);
    end

    %% H. Risk-aware controller

    feasibleA = pOutA(k)<=epsilon;
    feasibleB = pOutB(k)<=epsilon;

    bothInfeasible(k) = ~feasibleA && ~feasibleB;

    if feasibleA && feasibleB
        % Both meet the outage requirement:
        % select higher predicted mean rate.
        if meanRateA(k)>=meanRateB(k)
            choiceRisk(k)=1;
        else
            choiceRisk(k)=2;
        end

    elseif feasibleA
        choiceRisk(k)=1;

    elseif feasibleB
        choiceRisk(k)=2;

    else
        % No feasible single link exists.
        % Select the link with lower outage probability.
        if pOutA(k)<=pOutB(k)
            choiceRisk(k)=1;
        else
            choiceRisk(k)=2;
        end
    end

    if choiceRisk(k)==1
        rateRisk(k)=mean(RAt);
        outRisk(k)=mean(RAt<Rmin);
    else
        rateRisk(k)=mean(RBt);
        outRisk(k)=mean(RBt<Rmin);
    end

    %% I. Fresh-pilot-assisted decision

    % Optimistic benchmark:
    % Both links are measured at the future decision
    % instant. There is no processing/switching delay.
    %
    % The pilot observation is noisy.
    pilotNoiseA = sqrt(pilotNoiseVar/2)* ...
        (randn(Ntest,1)+1j*randn(Ntest,1));

    pilotNoiseB = sqrt(pilotNoiseVar/2)* ...
        (randn(Ntest,1)+1j*randn(Ntest,1));

    zA = hAt+pilotNoiseA;
    zB = hBt+pilotNoiseB;

    % Equal reference SNR: choose stronger measured power
    chooseA = abs(zA).^2 >= abs(zB).^2;

    pilotRateSamples = ...
        chooseA.*RAt + (~chooseA).*RBt;

    pilotOutSamples = ...
        chooseA.*(RAt<Rmin) + ...
        (~chooseA).*(RBt<Rmin);

    % Account for pilot overhead
    ratePilot(k) = ...
        (1-pilotOverhead)*mean(pilotRateSamples);

    outPilot(k) = mean(pilotOutSamples);

    %% J. Perfect-future-information benchmark

    rateOracle(k)=mean(max(RAt,RBt));

end

%% 9. Calculate expected decision regret

regretStale = rateOracle-rateStale;
regretPred = rateOracle-ratePred;
regretRisk = rateOracle-rateRisk;

% Additional value of a fresh pilot over prediction
pilotNetGain = ratePilot-ratePred;

%% 10. Figure 1: Physical channel and uncertainty

figure('Color','w','Position',[80 80 1200 750]);

subplot(2,2,1);

plot(timeMs,sampleA,'Color',[0.8 0.8 0.8]);
hold on;
plot(timeMs,sampleB,'Color',[0.75 0.85 0.95]);

plot(timeMs,meanPowerA,'r-','LineWidth',2);
plot(timeMs,meanPowerB,'b-','LineWidth',2);

grid on;
xlim([0 5]);
xlabel('Future time (ms)');
ylabel('Channel power');
title('Possible Future Multipath Channels');
legend('Link A samples','Link B samples', ...
    'Mean A','Mean B','Location','best');

subplot(2,2,2);

plot(timeMs,varChanA,'r-','LineWidth',2);
hold on;
plot(timeMs,varChanB,'b-','LineWidth',2);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Complex channel variance');
title('Future Channel Uncertainty');
legend('Link A','Link B','Location','best');

subplot(2,2,3);

plot(timeMs,meanRateA,'r-','LineWidth',2);
hold on;
plot(timeMs,meanRateB,'b-','LineWidth',2);
yline(Rmin,'k:','Required rate');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected rate (bit/s/Hz)');
title('Expected Future Spectral Efficiency');
legend('Link A','Link B','Required rate', ...
    'Location','best');

subplot(2,2,4);

plot(timeMs,pOutA,'r-','LineWidth',2);
hold on;
plot(timeMs,pOutB,'b-','LineWidth',2);
yline(epsilon,'k:','Outage limit');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Outage probability');
title('Predicted Reliability of Each Link');
legend('Link A','Link B','5% limit', ...
    'Location','best');

sgtitle('Experiment 4: Two Uncertain Multipath Links');

%% 11. Figure 2: Controller comparison

figure('Color','w','Position',[120 80 1200 850]);

subplot(3,2,1);

stairs(timeMs,choicePred,'b-','LineWidth',1.5);
hold on;
stairs(timeMs,choiceRisk,'m--','LineWidth',1.5);

grid on;
ylim([0.7 2.3]);
yticks([1 2]);
yticklabels({'Link A','Link B'});
xlabel('Prediction horizon (ms)');
title('Controller Link Selection');
legend('Predictive','Risk-aware','Location','best');

subplot(3,2,2);

plot(timeMs,outStale,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);
hold on;
plot(timeMs,outPred,'b-','LineWidth',1.5);
plot(timeMs,outRisk,'m--','LineWidth',1.5);
plot(timeMs,outPilot,'g-','LineWidth',1.5);
yline(epsilon,'r:','5% limit');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Actual outage probability');
title('Reliability of Selected Link');
legend('Stale','Predictive','Risk-aware', ...
    'Fresh pilot','Limit','Location','best');

subplot(3,2,3);

plot(timeMs,rateStale,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);
hold on;
plot(timeMs,ratePred,'b-','LineWidth',1.5);
plot(timeMs,rateRisk,'m--','LineWidth',1.5);
plot(timeMs,ratePilot,'g-','LineWidth',1.5);
plot(timeMs,rateOracle,'k:','LineWidth',2);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected rate (bit/s/Hz)');
title('Communication Performance');
legend('Stale','Predictive','Risk-aware', ...
    'Fresh pilot, net','Oracle','Location','best');

subplot(3,2,4);

plot(timeMs,regretStale,'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5);
hold on;
plot(timeMs,regretPred,'b-','LineWidth',1.5);
plot(timeMs,regretRisk,'m--','LineWidth',1.5);
yline(regretBudget,'r:','Regret budget');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected regret (bit/s/Hz)');
title('Decision Regret');
legend('Stale','Predictive','Risk-aware', ...
    'Regret budget','Location','best');

subplot(3,2,5);

plot(timeMs,pilotNetGain,'g-','LineWidth',2);
hold on;
yline(0,'k--');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Net rate gain (bit/s/Hz)');
title('Benefit of Fresh Channel Measurements');

subplot(3,2,6);

stairs(timeMs,double(bothInfeasible), ...
    'r-','LineWidth',1.5);

grid on;
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'Feasible link','Neither feasible'});
xlabel('Prediction horizon (ms)');
title('Can Either Link Meet the Reliability Target?');

sgtitle('Predictive and Risk-Aware Decisions');

%% 12. Decision lifetime from regret

regretMatrix = [regretStale regretPred regretRisk];
policyNames = {'Stale','Predictive','Risk-aware'};

fprintf('\n=== DECISION-REGRET LIFETIME ===\n');

for p = 1:3

    % Enforce budget throughout the prefix interval.
    crossing = find( ...
        cummax(regretMatrix(:,p))>regretBudget, ...
        1,'first');

    if isempty(crossing)
        lifetime = timeMs(end);
    elseif crossing==1
        lifetime = 0;
    else
        lifetime = timeMs(crossing-1);
    end

    fprintf('%-12s : %.3f ms\n', ...
        policyNames{p},lifetime);
end

%% 13. Summary results

fprintf('\n=== SYSTEM RELIABILITY ===\n');
fprintf('Reliability target: outage <= %.1f%%\n', ...
    100*epsilon);

fprintf('Neither link feasible at %.1f%% of horizons\n', ...
    100*mean(bothInfeasible));

fprintf('Fresh pilot beneficial at %.1f%% of horizons\n', ...
    100*mean(pilotNetGain>0));

fprintf('\n=== SELECTED HORIZONS ===\n');

checkMs = [0.5 1 2 5 10];

fprintf(['Horizon   Rate A   Rate B   ', ...
    'Out A   Out B   Risk choice\n']);

for j = 1:numel(checkMs)

    [~,idx] = min(abs(timeMs-checkMs(j)));

    fprintf('%5.1f ms   %5.2f    %5.2f   ', ...
        timeMs(idx),meanRateA(idx),meanRateB(idx));

    fprintf('%5.3f   %5.3f      Link %c\n', ...
        pOutA(idx),pOutB(idx), ...
        char('A'+choiceRisk(idx)-1));
end

fprintf('\nExperiment 4 complete.\n');
