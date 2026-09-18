function M = multivariateInformation(F,H)

%% =========================================================
% INPUT
%
% F = N x d
% H = N x d
%% =========================================================

F = double(F);
H = double(H);


%% =========================================================
% STEP 1
% GAUSSIAN COPULA TRANSFORMATION
%% =========================================================

Fg = gaussianize(F);

Hg = gaussianize(H);


%% =========================================================
% STEP 2
% COVARIANCE MATRICES
%% =========================================================

CF = cov(Fg);

CH = cov(Hg);

joint = [Fg Hg];

CJ = cov(joint);


%% =========================================================
% STEP 3
% SMALL REGULARIZATION
%% =========================================================

lambda = 1e-6;

CF = CF + lambda*eye(size(CF));

CH = CH + lambda*eye(size(CH));

CJ = CJ + lambda*eye(size(CJ));


%% =========================================================
% STEP 4
% LOG DETERMINANTS
%% =========================================================

logdetF = stableLogDet(CF);

logdetH = stableLogDet(CH);

logdetJ = stableLogDet(CJ);


%% =========================================================
% STEP 5
% MUTUAL INFORMATION
%
% I(F;H)
%% =========================================================

MI = ...
    0.5 * ...
    (logdetF + logdetH - logdetJ) ...
    / log(2);


MI = max(MI,0);


%% =========================================================
% STEP 6
% DIFFERENTIAL ENTROPY
%% =========================================================

dF = size(Fg,2);

dH = size(Hg,2);


HF = ...
    0.5 * ...
    ( ...
      dF*log(2*pi*exp(1)) + ...
      logdetF ...
    ) / log(2);


HH = ...
    0.5 * ...
    ( ...
      dH*log(2*pi*exp(1)) + ...
      logdetH ...
    ) / log(2);


%% =========================================================
% STEP 7
% NORMALIZED MUTUAL INFORMATION
%% =========================================================

NMI = ...
    2*MI/(HF + HH);


%% =========================================================
% STEP 8
% DIRECTIONAL SHARED INFORMATION
%% =========================================================

SharedFTIR = MI/HF;

SharedHE = MI/HH;


%% =========================================================
% KEEP NUMERIC RANGE SENSIBLE
%% =========================================================

NMI = max(0,min(1,NMI));

SharedFTIR = ...
    max(0,min(1,SharedFTIR));

SharedHE = ...
    max(0,min(1,SharedHE));


%% =========================================================
% OUTPUT
%% =========================================================

M.MI = MI;

M.HF = HF;

M.HH = HH;

M.NMI = NMI;

M.SharedFTIR = SharedFTIR;

M.SharedHE = SharedHE;

end