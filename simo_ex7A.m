
%% EXPERIMENT 7A: SPATIAL CHANNEL AGING

% Channel vector at time zero
hInitial = hA(:,1);

% Normalized similarity between initial and future channels
rhoSpatial = abs(hInitial'*hA) ./ ...
    (norm(hInitial)*sqrt(sum(abs(hA).^2,1)));

% Fraction of available array power captured by stale MRC
alignmentEfficiency = rhoSpatial.^2;

figure('Color','w');

plot(timeMs,rhoSpatial,'b-','LineWidth',2);
hold on;

plot(timeMs,alignmentEfficiency,'r--','LineWidth',2);

grid on;
xlabel('Prediction horizon (ms)');
ylabel('Spatial alignment');
title('How Quickly Does the SIMO Channel Direction Age?');

legend('Channel direction correlation', ...
       'Stale combining efficiency', ...
       'Location','best');

ylim([0 1.05]);

% Verify the matrix relationship
gainFromAlignment = ...
    sum(abs(hA).^2,1).*alignmentEfficiency;

fprintf('Maximum identity error: %.3e\n', ...
    max(abs(gainFromAlignment-gainAStale)));
