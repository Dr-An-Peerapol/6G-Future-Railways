
%% EXPERIMENT 2: PREDICTIVE EQUALIZATION WITH UNCERTAINTY
% Append after Experiment 1.
%
% Channel: y = h*x + n, with QPSK symbols.
%
% Compare:
%   1. Stale CSI equalizer
%   2. Predictive robust LMMSE equalizer
%   3. Predictive equalizer + pilot refresh
%   4. Perfect-CSI receiver (benchmark)
%
% Assumptions:
% - Path gains are known perfectly at t = 0.
% - Nominal Doppler frequencies are known.
% - Actual Dopplers differ by independent Gaussian errors.
% - Pilot refresh is instantaneous (optimistic benchmark).
% - Each horizon is an independent Monte Carlo experiment.

rng(21);

%% Parameters
sigmaF = 50;             % Doppler uncertainty std (Hz)
SNRrx_dB = 15;           % Receiver SNR
N0 = 10^(-SNRrx_dB/10);

Nmc = 5000;              % Monte Carlo trials per horizon
horizonMs = (0:0.1:10).';
Nh = numel(horizonMs);

NpilotAvg = 8;           % Pilot samples per refresh
refreshThreshold = 0.15;

% Result arrays: rows = horizons, columns = channels
predVar = zeros(Nh,2);
empVar = zeros(Nh,2);
refresh = false(Nh,2);

berStale = zeros(Nh,2);
berPred  = zeros(Nh,2);
berPilot = zeros(Nh,2);
berIdeal = zeros(Nh,2);

mseNominal = zeros(Nh,2);
mseRobust  = zeros(Nh,2);

%% Simulate both channel types
for ch = 1:2

    fNom = fList{ch}(:).';
    a = aList{ch}(:);
    L = numel(fNom);

    % Current channel is known at time t = 0.
    h0 = sum(a);

    for ii = 1:Nh

        tau = horizonMs(ii)*1e-3;

        %% A. Analytic future-channel uncertainty

        % Attenuation of the conditional complex mean
        coherenceFactor = ...
            exp(-0.5*(2*pi*sigmaF*tau)^2);

        % Nominal future channel
        hNom = exp(1j*2*pi*tau*fNom)*a;

        % Conditional mean and error variance
        mu = coherenceFactor*hNom;

        vh = (1 - coherenceFactor^2)*sum(abs(a).^2);

        predVar(ii,ch) = vh;

        %% B. Monte Carlo true future channels

        % Different possible Doppler values
        fActual = repmat(fNom,Nmc,1) ...
                + sigmaF*randn(Nmc,L);

        % True channel for each realization
        hActual = exp(1j*2*pi*tau*fActual)*a;

        % Verify the analytical uncertainty prediction
        empVar(ii,ch) = mean(abs(hActual-mu).^2);

        %% C. Transmit random QPSK data

        bitI = 2*randi([0 1],Nmc,1)-1;
        bitQ = 2*randi([0 1],Nmc,1)-1;

        x = (bitI+1j*bitQ)/sqrt(2);

        noise = sqrt(N0/2)* ...
            (randn(Nmc,1)+1j*randn(Nmc,1));

        y = hActual.*x+noise;

        %% D. Receiver 1: stale CSI

        wStale = conj(h0)/(abs(h0)^2+N0);
        zStale = wStale*y;

        %% E. Receiver 2: predictive robust LMMSE

        wRobust = conj(mu)/(abs(mu)^2+vh+N0);
        zPred = wRobust*y;

        %% F. Conventional predicted-channel MMSE

        % This receiver ignores channel prediction uncertainty.
        wNominal = conj(hNom)/(abs(hNom)^2+N0);
        zNominal = wNominal*y;

        % Compare equalized-symbol MSE
        mseNominal(ii,ch) = mean(abs(zNominal-x).^2);
        mseRobust(ii,ch)  = mean(abs(zPred-x).^2);

        %% G. Decide whether a new pilot is needed

        % A heuristic uncertainty indicator.
        uncertaintyRatio = ...
            vh/(abs(mu)^2+vh+N0);

        needPilot = uncertaintyRatio > refreshThreshold;
        refresh(ii,ch) = needPilot;

        if needPilot
            % Fresh pilot measurement at the target time.
            % Average NpilotAvg independent observations,
            % assuming the channel remains constant during them.
            pilotNoise = sqrt(N0/(2*NpilotAvg))* ...
                (randn(Nmc,1)+1j*randn(Nmc,1));

            hPilot = hActual+pilotNoise;

            wPilot = conj(hPilot)./(abs(hPilot).^2+N0);

            zPilot = wPilot.*y;
        else
            zPilot = zPred;
        end

        %% H. Perfect CSI benchmark

        wIdeal = conj(hActual)./(abs(hActual).^2+N0);
        zIdeal = wIdeal.*y;

        %% I. QPSK bit error rates

        % Compare the real and imaginary bit decisions.
        countBER = @(z) ...
            (nnz((real(z)>=0) ~= (real(x)>=0)) + ...
             nnz((imag(z)>=0) ~= (imag(x)>=0))) ...
             /(2*Nmc);

        berStale(ii,ch) = countBER(zStale);
        berPred(ii,ch)  = countBER(zPred);
        berPilot(ii,ch) = countBER(zPilot);
        berIdeal(ii,ch) = countBER(zIdeal);

    end
end

%% Figure 1: BER and uncertainty
figure('Color','w','Position',[100 100 1200 750]);

for ch = 1:2

    subplot(2,2,ch);

    % Floor only for logarithmic plotting of zero observed errors
    berFloor = 1/(2*Nmc);

    semilogy(horizonMs,max(berStale(:,ch),berFloor), ...
        'Color',[0.5 0.5 0.5],'LineWidth',1.5);
    hold on;

    semilogy(horizonMs,max(berPred(:,ch),berFloor), ...
        'b-','LineWidth',1.5);

    semilogy(horizonMs,max(berPilot(:,ch),berFloor), ...
        'm-','LineWidth',1.5);

    semilogy(horizonMs,max(berIdeal(:,ch),berFloor), ...
        'k--','LineWidth',1.5);

    grid on;
    xlabel('Prediction horizon (ms)');
    ylabel('QPSK bit error rate');
    title([names{ch}, ': BER']);
    legend('Stale CSI','Predictive LMMSE', ...
           'Predict + Pilot','Perfect CSI', ...
           'Location','best');

    subplot(2,2,ch+2);

    plot(horizonMs,predVar(:,ch), ...
        'b-','LineWidth',2);
    hold on;

    plot(horizonMs,empVar(:,ch), ...
        'k:','LineWidth',1.5);

    % Mark horizons at which the receiver requests pilots
    idx = find(refresh(:,ch));

    scatter(horizonMs(idx),predVar(idx,ch), ...
        15,'m','filled');

    grid on;
    xlabel('Prediction horizon (ms)');
    ylabel('Future channel error variance');
    title([names{ch}, ': Prediction Uncertainty']);

    legend('Analytical variance','Monte Carlo variance', ...
        'Pilot requested','Location','best');
end

sgtitle('Uncertainty-Aware Predictive Equalization');

%% Figure 2: Why robust LMMSE differs from nominal MMSE

figure('Color','w','Position',[150 150 1000 350]);

for ch = 1:2

    subplot(1,2,ch);

    plot(horizonMs,10*log10(mseNominal(:,ch)), ...
        'r--','LineWidth',1.5);
    hold on;

    plot(horizonMs,10*log10(mseRobust(:,ch)), ...
        'b-','LineWidth',1.5);

    grid on;
    xlabel('Prediction horizon (ms)');
    ylabel('Equalized-symbol MSE (dB)');
    title(names{ch});

    legend('Nominal MMSE','Uncertainty-aware LMMSE', ...
        'Location','best');
end

sgtitle('Effect of Channel Uncertainty on Equalizer MSE');

%% Selected numerical results

fprintf('\n=== EXPERIMENT 2 RESULTS ===\n');

checkMs = [1 5 10];

for ch = 1:2

    fprintf('\n%s\n',names{ch});

    fprintf(['Horizon   Var(h)   BER stale   ', ...
             'BER predict   BER + pilot   BER ideal\n']);

    for k = 1:numel(checkMs)

        [~,idx] = min(abs(horizonMs-checkMs(k)));

        fprintf('%5.1f ms   %.4f   %.5f      %.5f', ...
            horizonMs(idx),predVar(idx,ch), ...
            berStale(idx,ch),berPred(idx,ch));

        fprintf('        %.5f      %.5f\n', ...
            berPilot(idx,ch),berIdeal(idx,ch));
    end

    fprintf('Pilot-request fraction across horizons: %.1f%%\n', ...
        100*mean(refresh(:,ch)));
end
