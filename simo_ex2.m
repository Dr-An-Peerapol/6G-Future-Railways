
%% SIMO Lesson 2: Coherent Spatial Combining

% Transmitted QPSK symbol
x = (1 + 1j)/sqrt(2);

% Noise-free signals received by four antennas
y = a*x;

% Naive addition
naiveOutput = sum(y);

% Coherent combining using the steering vector
coherentOutput = (a'*y)/(a'*a);

fprintf('Transmitted symbol:       %.4f %+.4fj\n', ...
    real(x), imag(x));

fprintf('Naively combined output:  %.4f %+.4fj\n', ...
    real(naiveOutput), imag(naiveOutput));

fprintf('Coherently recovered:     %.4f %+.4fj\n', ...
    real(coherentOutput), imag(coherentOutput));
