function X_processed = preprocess_FTIR(X, polyOrder, frameLength)

% =========================================================
% Input:
% X = N pixels × 965 spectral channels
%
% Processing:
% 1. Shift minimum of each spectrum to zero
% 2. L2 normalization
% 3. Savitzky-Golay second derivative
% =========================================================


%% ---------------------------------------------------------
% 1. Shift every spectrum so minimum = 0
%% ---------------------------------------------------------

minimumValue = min(X, [], 2);

X_shift = X - minimumValue;


%% ---------------------------------------------------------
% 2. Normalize every spectrum to length = 1
%% ---------------------------------------------------------

spectrumNorm = sqrt(sum(X_shift.^2, 2));

% Prevent division by zero
spectrumNorm(spectrumNorm == 0) = 1;

X_norm = X_shift ./ spectrumNorm;


%% ---------------------------------------------------------
% 3. Savitzky-Golay second derivative
%
% Polynomial order = 5
% Window = 17
%% ---------------------------------------------------------

[~, G] = sgolay(polyOrder, frameLength);

derivativeOrder = 2;

SGcoeff = factorial(derivativeOrder) * ...
          G(:, derivativeOrder + 1);


% Apply along spectral dimension
X_processed = conv2(X_norm, SGcoeff', 'same');


%% ---------------------------------------------------------
% Remove unreliable SG boundary regions
%% ---------------------------------------------------------

halfWindow = (frameLength - 1) / 2;

X_processed(:, 1:halfWindow) = NaN;

X_processed(:, end-halfWindow+1:end) = NaN;

end