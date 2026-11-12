%% EXPERIMENT 1: SINGLE-PATH vs TWO-PATH CHANNEL PREDICTABILITY
% Research objective:
% Investigate how Doppler-model errors affect future channel prediction.
%
% Scenario: 28 GHz carrier, train moving at 350 km/h.
% No RIS and no EKF yet.
%
% Three predictors:
% 1. Hold the most recent pilot observation
% 2. Predict using correct Doppler frequencies
% 3. Predict using slightly incorrect Doppler frequencies
%
% All predictors use the same noisy pilot measurements.

clear; close all; clc;
rng(7);  % Reproducible simulation

%% 1. Physical parameters
c = 299792458;        % Speed of light (m/s)
fc = 28e9;            % Carrier frequency (Hz)
v = 350 / 3.6;        % Train speed (m/s)
lambda = c / fc;

fDmax = v / lambda;   % Maximum Doppler magnitude

fprintf('Train speed       : %.3f m/s\n', v);
fprintf('Wavelength        : %.3f mm\n', lambda * 1e3);
fprintf('Max Doppler       : %.2f Hz\n\n', fDmax);

%% 2. Simulation settings
SNRdB = 20;

Tobs = 2e-3;          % 2 ms of past observations
Tpilot = 20e-6;       % Pilot every 20 microseconds

Tfuture = 10e-3;      % Predict 10 ms into the future
dt = 5e-6;            % Fine time grid for plotting

tObs = (-Tobs:Tpilot:0).';
tFuture = (0:dt:Tfuture).';

% The prediction starts at t = 0.
% Past measurements occur at t <= 0.
% No future measurements are used.

%% 3. Define two physical channels

% Channel 1: one dominant propagation path
fList{1} = fDmax;
aList{1} = 1;

% Channel 2: two propagation paths
% Their Doppler difference is 1000 Hz.
fList{2} = [fDmax, fDmax - 1000];

% Equal path amplitudes, normalized average power
aList{2} = [1; 1] / sqrt(2);

names = {'Single-path channel', 'Two-path channel'};

% Deliberate Doppler model error
% Single path: +50 Hz
% Two paths: +50 Hz and -50 Hz
fBias{1} = 50;
fBias{2} = [50, -50];

%% 4. Estimate path gains and predict future channels

for k = 1:2

    fTrue = fList{k};
    aTrue = aList{k};

    % True channel at pilot observation times
    Aobs = exp(1j * 2*pi * (tObs * fTrue));
    hObs = Aobs * aTrue;

    % Complex Gaussian pilot noise
    signalPower = mean(abs(hObs).^2);
    noisePower = signalPower / (10^(SNRdB/10));

    noise = sqrt(noisePower/2) * ...
        (randn(size(hObs)) + 1j*randn(size(hObs)));

    yObs = hObs + noise;

    % Ground-truth future channel
    Afuture = exp(1j * 2*pi * (tFuture * fTrue));
    hTrue = Afuture * aTrue;

    % ---- Predictor 1: Last-observation hold ----
    hHold = yObs(end) * ones(size(tFuture));

    % ---- Predictor 2: Accurate Doppler model ----
    % Frequencies are supplied correctly.
    % Complex path gains are estimated by LS.
    aEst = Aobs \ yObs;
    hAccurate = Afuture * aEst;

    % ---- Predictor 3: Inaccurate Doppler model ----
    % The assumed frequencies contain small errors.
    fWrong = fTrue + fBias{k};

    AobsWrong = exp(1j * 2*pi * (tObs * fWrong));
    aEstWrong = AobsWrong \ yObs;

    AfutureWrong = exp(1j * 2*pi * ...
        (tFuture * fWrong));

    hWrong = AfutureWrong * aEstWrong;

    % Save all results
    results(k).true = hTrue;
    results(k).hold = hHold;
    results(k).accurate = hAccurate;
    results(k).wrong = hWrong;

    % Cumulative normalized mean-squared error:
    % NMSE(T) = sum_{t<=T}|h_pred-h_true|^2 /
    %           sum_{t<=T}|h_true|^2

    denominator = max(cumsum(abs(hTrue).^2), eps);

    results(k).nmseHold = ...
        cumsum(abs(hHold - hTrue).^2) ./ denominator;

    results(k).nmseAccurate = ...
        cumsum(abs(hAccurate - hTrue).^2) ./ denominator;

    results(k).nmseWrong = ...
        cumsum(abs(hWrong - hTrue).^2) ./ denominator;

end

%% 5. Plot received power and prediction error

figure('Color','w','Position',[100 100 1200 700]);

timeMs = tFuture * 1e3;

for k = 1:2

    % Upper row: received channel power
    subplot(2,2,k);

    plot(timeMs, abs(results(k).true).^2, ...
        'k-', 'LineWidth', 2); hold on;

    plot(timeMs, abs(results(k).accurate).^2, ...
        'b--', 'LineWidth', 1.5);

    plot(timeMs, abs(results(k).wrong).^2, ...
        'r-.', 'LineWidth', 1.5);

    plot(timeMs, abs(results(k).hold).^2, ...
        'Color',[0.55 0.55 0.55], ...
        'LineStyle',':', 'LineWidth',1.5);

    grid on;
    xlim([0 2]);
    xlabel('Future time (ms)');
    ylabel('|h(t)|^2');
    title([names{k}, ': Received Power']);

    legend('True channel','Correct Doppler', ...
        'Wrong Doppler','Last observation', ...
        'Location','best');

    % Lower row: cumulative channel prediction NMSE
    subplot(2,2,k+2);

    plot(timeMs, 10*log10(results(k).nmseHold + eps), ...
        'Color',[0.55 0.55 0.55], ...
        'LineWidth',1.5); hold on;

    plot(timeMs, 10*log10(results(k).nmseAccurate + eps), ...
        'b-', 'LineWidth',1.5);

    plot(timeMs, 10*log10(results(k).nmseWrong + eps), ...
        'r-', 'LineWidth',1.5);

    grid on;
    xlabel('Prediction horizon (ms)');
    ylabel('Cumulative NMSE (dB)');
    title([names{k}, ': Prediction Error']);

    legend('Last observation','Correct Doppler', ...
        'Wrong Doppler','Location','best');

end

sgtitle('Channel Predictability at 350 km/h, 28 GHz');

%% 6. Numerical results at selected prediction horizons

horizonsMs = [1 5 10];

fprintf('\n=== CUMULATIVE PREDICTION NMSE (dB) ===\n');

for k = 1:2

    fprintf('\n%s\n', names{k});
    fprintf('Horizon     Hold      Correct fD     Wrong fD\n');

    for m = 1:length(horizonsMs)

        idx = round(horizonsMs(m)*1e-3/dt) + 1;

        e1 = 10*log10(results(k).nmseHold(idx) + eps);
        e2 = 10*log10(results(k).nmseAccurate(idx) + eps);
        e3 = 10*log10(results(k).nmseWrong(idx) + eps);

        fprintf('%5.1f ms   %8.2f   %10.2f   %10.2f\n', ...
            horizonsMs(m), e1, e2, e3);
    end
end

%% 7. Theoretical two-path fading check

deltaFD = abs(diff(fList{2}));
firstNull = 1 / (2*deltaFD);

fprintf('\nTwo-path Doppler separation: %.0f Hz\n', deltaFD);
fprintf('First ideal fading null    : %.3f ms\n', ...
    firstNull*1e3);
