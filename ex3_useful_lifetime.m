
%% EXPERIMENT 3
% USEFUL LIFETIME AND VALUE OF CHANNEL INFORMATION
%
% Research questions:
% 1. How quickly does future channel uncertainty grow?
% 2. How long does old information remain useful for decisions?
% 3. Can prediction extend decision-relevant information lifetime?
% 4. When is acquiring another pilot worth its overhead?
%
% No RIS, EKF, or equalizer in this experiment.
% Each horizon is evaluated as a separate future decision.
%
% MATLAB: base functions only.

clear; close all; clc;
rng(14);

%% 1. Physical parameters

c = 299792458;
fc = 28e9;
v = 350/3.6;
lambda = c/fc;
fDmax = v/lambda;

fprintf('Train speed: %.2f km/h\n',v*3.6);
fprintf('Maximum Doppler: %.2f Hz\n',fDmax);

% Two paths arriving at different angles
% Nominal Dopplers stay below maximum radial Doppler.
fNom = [0.9*fDmax, 0.9*fDmax-1000];

sigmaFD = 50;           % Doppler error std, Hz

% Link A: equal-amplitude two-path channel
a = [1;1]/sqrt(2);

% Link B: constant-power single-path alternative
PB = 0.8;

%% 2. Communication and research settings

snrLinear = 10;         % 10 dB reference SNR
Rmin = 2.5;             % Minimum reliable rate, bit/s/Hz
epsilon = 0.05;         % Maximum permitted outage probability

% Maximum allowed expected regret
regretBudget = 0.10;    % bit/s/Hz

% Hypothetical future pilot
pilotNoiseVar = 0.05;   % Complex observation noise variance
pilotOverhead = 0.05;   % 5 percent of a transmission slot

% Prediction horizons
timeMs = (0:0.05:10).';
timeSec = timeMs*1e-3;
Nt = length(timeMs);

% Monte Carlo settings
Nmc = 2500;
Nprobe = 2000;


% Link B has a known, constant power.
RB = log2(1+snrLinear*PB);

%% 3. Draw possible physical channel trajectories
% An individual Doppler error remains constant over 10 ms.
% The errors of the two paths are independent.
%
% Crucially, the system does not know their actual values.
% It knows only the nominal values and sigmaFD.

fActual = repmat(fNom,Nmc,1) ...
        + sigmaFD*randn(Nmc,2);

% Independent samples for hypothetical pilot observations
fProbe = repmat(fNom,Nprobe,1) ...
       + sigmaFD*randn(Nprobe,2);

%% 4. Allocate arrays

meanPowerA = zeros(Nt,1);
powerVariance = zeros(Nt,1);

channelErrorVar = zeros(Nt,1);
channelErrorVarMC = zeros(Nt,1);
complexNMSE = zeros(Nt,1);

expectedRateA = zeros(Nt,1);
outageA = zeros(Nt,1);

rateStale = zeros(Nt,1);
ratePredict = zeros(Nt,1);
rateRisk = zeros(Nt,1);
rateOracle = zeros(Nt,1);

outStale = zeros(Nt,1);
outPredict = zeros(Nt,1);
outRisk = zeros(Nt,1);

choicePredict = zeros(Nt,1);
choiceRisk = zeros(Nt,1);

EVSI = zeros(Nt,1);
EVPI = zeros(Nt,1);
netPilotGain = zeros(Nt,1);
rateWithPilotOption = zeros(Nt,1);

% Store a few sample trajectories for plotting
Nshow = 12;
samplePowers = zeros(Nshow,Nt);

%% 5. Main simulation loop

for k = 1:Nt

    tau = timeSec(k);

    %% A. True future channels (Monte Carlo)

    hA = ( ...
        exp(1j*2*pi*tau*fActual(:,1)) + ...
        exp(1j*2*pi*tau*fActual(:,2)) ...
        )/sqrt(2);

    PA = abs(hA).^2;
    samplePowers(:,k) = PA(1:Nshow);

    RA = log2(1+snrLinear*PA);

    % Link B rate remains constant.
    avgRA = mean(RA);

    expectedRateA(k) = avgRA;
    outageA(k) = mean(RA < Rmin);

    %% B. Analytic future channel distribution moments

    % Conditional mean: E[h_A(tau) | current information]
    coherence = exp(-0.5*(2*pi*sigmaFD*tau)^2);

    hNom = ( ...
        exp(1j*2*pi*tau*fNom(1)) + ...
        exp(1j*2*pi*tau*fNom(2)) ...
        )/sqrt(2);

    muH = coherence*hNom;

    % Exact channel prediction error variance
    % for independent Gaussian Doppler errors
    varH = (1-coherence^2)*sum(abs(a).^2);

    channelErrorVar(k) = varH;

    % Monte Carlo check
    channelErrorVarMC(k) = mean(abs(hA-muH).^2);

    %% C. Mean received power

    deltaFD = fNom(1)-fNom(2);

    relativeCoherence = exp( ...
        -0.5*(2*pi*tau)^2*(2*sigmaFD^2));

    meanPA = 1 + relativeCoherence * ...
        cos(2*pi*deltaFD*tau);

    meanPowerA(k) = meanPA;

    powerVariance(k) = mean((PA-meanPA).^2);

    complexNMSE(k) = ...
        mean(abs(hA-muH).^2)/mean(PA);

    %% D. Strategy 1: Stale channel information

    % At t = 0, Link A has power 2, so it is selected.
    % The stale policy never changes this initial choice.
    rateStale(k) = avgRA;
    outStale(k) = outageA(k);

    %% E. Strategy 2: Predictive decision

    if avgRA > RB
        choicePredict(k) = 1;  % Choose Link A
        ratePredict(k) = avgRA;
        outPredict(k) = outageA(k);
    else
        choicePredict(k) = 0;  % Choose Link B
        ratePredict(k) = RB;
        outPredict(k) = 0;
    end

    %% F. Strategy 3: Risk-aware decision

    % Choose A only when:
    % 1) Its expected rate exceeds that of B.
    % 2) Its outage probability is below epsilon.
    %
    % B is assumed to satisfy the reliability constraint.

    if (avgRA > RB) && (outageA(k) <= epsilon)
        choiceRisk(k) = 1;
        rateRisk(k) = avgRA;
        outRisk(k) = outageA(k);
    else
        choiceRisk(k) = 0;
        rateRisk(k) = RB;
        outRisk(k) = 0;
    end

    %% G. Perfect future information benchmark

    % Oracle knows actual Link A power before choosing.
    % Not a realizable policy without measurement cost.
    oracleSamples = max(RA,RB);
    rateOracle(k) = mean(oracleSamples);

    % Value of perfect information (upper benchmark)
    currentBest = max(avgRA,RB);
    EVPI(k) = rateOracle(k)-currentBest;

    %% H. Hypothetical new pilot measurement

    % Generate future Link A states.
    hProbe = ( ...
        exp(1j*2*pi*tau*fProbe(:,1)) + ...
        exp(1j*2*pi*tau*fProbe(:,2)) ...
        )/sqrt(2);

    % Observe z = h_A + complex Gaussian noise.
    pilotNoise = sqrt(pilotNoiseVar/2)* ...
        (randn(Nprobe,1)+1j*randn(Nprobe,1));

    z = hProbe+pilotNoise;

    %% I. Bayesian posterior approximation

    % Each Monte Carlo hA sample is a prior particle.
    %
    % Likelihood:
    % p(z|hA) proportional to
    % exp(-abs(z-hA)^2/pilotNoiseVar)
    %
    % Matrix size: Nmc x Nprobe

    residual = bsxfun(@minus,hA,z.');

    logWeights = ...
        -(abs(residual).^2)/pilotNoiseVar;

    % Numerically stable exponential normalization
    logWeights = bsxfun( ...
        @minus,logWeights,max(logWeights,[],1));

    weights = exp(logWeights);

    weights = bsxfun( ...
        @rdivide,weights,sum(weights,1));

    % Expected rate for Link A after observing each z
    postMeanRA = RA.'*weights;

    % After receiving z, choose A or B.
    afterPilotRate = mean(max(postMeanRA,RB));

    % Expected oracle rate, evaluated using the SAME
    % posterior weights and prospective observations
    oraclePosterior = max(RA, RB).' * weights;
    oracleRateAligned = mean(oraclePosterior);

    % Gross expected value of sample information
    % Small negative MC fluctuations are clipped to zero.
    EVSI(k) = max(afterPilotRate-currentBest,0);
    EVPI(k) = max(0, oracleRateAligned - currentBest);
    % Charge a fraction of the slot for the pilot.
    netPilotGain(k) = ...
        (1-pilotOverhead)*afterPilotRate-currentBest;

    % Best expected rate when pilot use is optional
    rateWithPilotOption(k) = max( ...
        currentBest, ...
        (1-pilotOverhead)*afterPilotRate);

end

%% 6. Decision regret and useful lifetime

regretStale = rateOracle-rateStale;
regretPredict = rateOracle-ratePredict;
regretRisk = rateOracle-rateRisk;

% Lifetime = last grid point before expected regret
% first exceeds the allowed budget.
%
% Prefix requirement makes lifetime meaningful even
% when the fading behavior is nonmonotonic.

idxStale = find( ...
    cummax(regretStale)>regretBudget,1,'first');

idxPredict = find( ...
    cummax(regretPredict)>regretBudget,1,'first');

if isempty(idxStale)
    lifetimeStale = timeMs(end);
else
    lifetimeStale = timeMs(max(1,idxStale-1));
end

if isempty(idxPredict)
    lifetimePredict = timeMs(end);
else
    lifetimePredict = timeMs(max(1,idxPredict-1));
end

%% 7. Figure 1: Channel and uncertainty

figure('Color','w','Position',[80 80 1250 450]);

subplot(1,3,1);

plot(timeMs,samplePowers, ...
    'Color',[0.78 0.78 0.78]); hold on;

plot(timeMs,meanPowerA,'b-', ...
    'LineWidth',2);

yline(PB,'k--','Link B');

grid on;
xlim([0 2]);
xlabel('Prediction horizon (ms)');
ylabel('Received channel power');
title('Possible Future Channels');

subplot(1,3,2);

plot(timeMs,channelErrorVar,'b-', ...
    'LineWidth',2); hold on;

plot(timeMs,channelErrorVarMC,'k:', ...
    'LineWidth',1.5);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Complex channel error variance');
title('Uncertainty Calibration');

legend('Analytical','Monte Carlo', ...
    'Location','southeast');

subplot(1,3,3);

plot(timeMs, ...
    10*log10(max(complexNMSE,1e-7)), ...
    'r-','LineWidth',2);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Complex prediction NMSE (dB)');
title('Aging of Complex CSI');

sgtitle('Channel Information Under Doppler Uncertainty');

%% 8. Figure 2: Decisions and reliability

figure('Color','w','Position',[100 80 1250 850]);

subplot(3,2,1);

stairs(timeMs,choicePredict, ...
    'b-','LineWidth',1.5); hold on;

stairs(timeMs,choiceRisk, ...
    'm--','LineWidth',1.5);

grid on;
xlim([0 2]);
ylim([-0.2 1.2]);

yticks([0 1]);
yticklabels({'Link B','Link A'});

xlabel('Prediction horizon (ms)');
title('Chosen Communication Link');

legend('Predictive','Risk-aware', ...
    'Location','best');

subplot(3,2,2);

plot(timeMs,outStale, ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5); hold on;

plot(timeMs,outPredict, ...
    'b-','LineWidth',1.5);

plot(timeMs,outRisk, ...
    'm--','LineWidth',1.5);

yline(epsilon,'r:','Outage limit');

grid on;
ylim([0 1]);
xlim([0 2]);
xlabel('Prediction horizon (ms)');
ylabel('Outage probability');
title('Communication Reliability');

legend('Stale','Predictive','Risk-aware', ...
    'Outage limit','Location','best');

subplot(3,2,3);

plot(timeMs,ratePredict, ...
    'b-','LineWidth',1.5); hold on;

plot(timeMs,rateWithPilotOption, ...
    'g-','LineWidth',1.5);

plot(timeMs,rateOracle, ...
    'k--','LineWidth',1.5);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected rate (bit/s/Hz)');
title('Does New Information Improve Decisions?');

legend('Predictive','Pilot optional','Oracle', ...
    'Location','best');

subplot(3,2,4);

plot(timeMs,regretStale, ...
    'Color',[0.5 0.5 0.5], ...
    'LineWidth',1.5); hold on;

plot(timeMs,regretPredict, ...
    'b-','LineWidth',1.5);

plot(timeMs,regretRisk, ...
    'm--','LineWidth',1.5);

yline(regretBudget,'r:','Regret budget');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected regret (bit/s/Hz)');
title('Decision-Relevant Information Lifetime');

legend('Stale','Predictive','Risk-aware', ...
    'Regret budget','Location','best');

subplot(3,2,5);

plot(timeMs,EVSI,'b-','LineWidth',1.5);
hold on;

plot(timeMs,EVPI,'k--','LineWidth',1.5);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Expected rate gain (bit/s/Hz)');
title('The Value of Information');

legend('Noisy future pilot','Perfect information', ...
    'Location','best');

subplot(3,2,6);

plot(timeMs,netPilotGain, ...
    'g-','LineWidth',1.5);
hold on;

yline(0,'k--');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Net gain (bit/s/Hz)');
title('When Is Another Pilot Worthwhile?');

sgtitle('Information Lifetime and Decision Reliability');

figure('Color','w');

plot(timeMs, ratePredict - RB, ...
    'b-', 'LineWidth', 2);

hold on;
yline(0, 'k--');

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Rate gain over fixed Link B (bit/s/Hz)');
title('Decision Benefit of Predictive Link Selection');

%% 9. Print research results

fprintf('\n=== DECISION-RELEVANT INFORMATION LIFETIME ===\n');
fprintf('Allowed regret: %.3f bit/s/Hz\n',regretBudget);
fprintf('Stale policy     : %.2f ms\n',lifetimeStale);
fprintf('Predictive policy: %.2f ms\n',lifetimePredict);

fprintf('\n=== SELECTED PREDICTION HORIZONS ===\n');

checkTimes = [0.5 1 2 5 10];

fprintf(['Horizon    CSI-NMSE    Pout(A)   ', ...
    'Regret(stale)   Regret(pred)   Pilot net\n']);

for i = 1:length(checkTimes)

    [~,idx] = min(abs(timeMs-checkTimes(i)));

    fprintf('%5.1f ms    %7.2f dB    %6.3f', ...
        timeMs(idx), ...
        10*log10(max(complexNMSE(idx),1e-12)), ...
        outageA(idx));

    fprintf('      %8.3f      %8.3f      %+8.3f\n', ...
        regretStale(idx), ...
        regretPredict(idx), ...
        netPilotGain(idx));
end

fprintf('\nPilot worthwhile at %.1f%% of tested horizons.\n', ...
    100*mean(netPilotGain>0));

fprintf('\nExperiment complete.\n');
