function Xout = preprocessFTIR(X)

%% =========================================================
% FTIR SPECTRAL PREPROCESSING
%
% Input:
% X = N pixels x 965 spectral channels
%
% Output:
% Xout = N x 965
%
% Processing:
% 1. Shift minimum to zero
% 2. L2 normalization
% 3. Savitzky-Golay second derivative
%
% SG parameters:
% Polynomial order = 5
% Window length    = 17
%% =========================================================


X = double(X);


%% =========================================================
% STEP 1
% SHIFT EACH SPECTRUM SO MINIMUM = 0
%% =========================================================

minimumValue = min(X,[],2);

X = X - minimumValue;


%% =========================================================
% STEP 2
% L2 NORMALIZATION
%% =========================================================

normValue = sqrt(sum(X.^2,2));

% Prevent division by zero
normValue(normValue < eps) = 1;

X = X ./ normValue;


%% =========================================================
% STEP 3
% DEFINE CONTINUOUS SPECTRAL REGIONS
%
% Original retained regions were:
%
% 13:203
% 282:675
% 780:908
% 1065:1315
%
% After reducing to 965 channels:
%
% Region 1 = 1:191
% Region 2 = 192:585
% Region 3 = 586:714
% Region 4 = 715:965
%% =========================================================

regions = {
    1:191
    192:585
    586:714
    715:965
    };


%% =========================================================
% STEP 4
% SAVITZKY-GOLAY PARAMETERS
%% =========================================================

polyOrder = 5;

frameLength = 17;

derivativeOrder = 2;


%% =========================================================
% CREATE SAVITZKY-GOLAY DIFFERENTIATION FILTER
%% =========================================================

[~,g] = sgolay(polyOrder,frameLength);

% MATLAB documentation:
%
% derivative filter =
%
% factorial(p) / (-dt)^p * g(:,p+1)
%
% Here spectral sample spacing is treated as 1 index.
%
% For second derivative:
% (-1)^2 = 1
%% =========================================================

dt = 1;

SGkernel = ...
    factorial(derivativeOrder) / ...
    ((-dt)^derivativeOrder) * ...
    g(:,derivativeOrder+1);


%% =========================================================
% OUTPUT MATRIX
%% =========================================================

Xout = zeros(size(X));


%% =========================================================
% STEP 5
% APPLY SECOND DERIVATIVE TO EACH CONTINUOUS REGION
%% =========================================================

halfWindow = floor(frameLength/2);


for r = 1:length(regions)

    idx = regions{r};

    regionData = X(:,idx);


    %% -----------------------------------------------------
    % Apply SG second derivative along spectral dimension
    %
    % regionData:
    %
    % N pixels x spectral channels
    %
    % SGkernel' operates along columns.
    %% -----------------------------------------------------

    derivativeData = ...
        conv2(regionData,SGkernel','same');


    %% -----------------------------------------------------
    % EDGE HANDLING
    %
    % Convolution near the first/last 8 spectral channels
    % of each continuous region is unreliable because the
    % full 17-point window is unavailable.
    %
    % Replace those edge values by the nearest reliable
    % derivative value.
    %% -----------------------------------------------------

    derivativeData(:,1:halfWindow) = ...
        repmat( ...
        derivativeData(:,halfWindow+1), ...
        1,halfWindow);


    derivativeData(:,end-halfWindow+1:end) = ...
        repmat( ...
        derivativeData(:,end-halfWindow), ...
        1,halfWindow);


    %% -----------------------------------------------------
    % Store
    %% -----------------------------------------------------

    Xout(:,idx) = derivativeData;

end


%% =========================================================
% REMOVE POSSIBLE NUMERICAL NON-FINITE VALUES
%% =========================================================

Xout(~isfinite(Xout)) = 0;

end