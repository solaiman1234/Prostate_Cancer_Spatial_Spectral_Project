function metrics = ...
    calculateInformationMetrics(F,H,numBins)

%% =========================================================
% F = FTIR energy values
% H = H&E grayscale values
%
% Both should already be normalized to [0,1]
%% =========================================================


%% =========================================================
% Histogram edges
%% =========================================================

edges = linspace(0,1,numBins+1);


%% =========================================================
% Joint histogram
%% =========================================================

jointCounts = histcounts2( ...
    F,H,edges,edges);


totalCount = sum(jointCounts(:));


if totalCount == 0

    error('No valid samples for MI calculation.');

end


Pjoint = jointCounts / totalCount;


%% =========================================================
% Marginal probabilities
%% =========================================================

PF = sum(Pjoint,2);

PH = sum(Pjoint,1);


%% =========================================================
% FTIR entropy
%% =========================================================

p = PF(PF > 0);

HF = -sum(p .* log2(p));


%% =========================================================
% H&E entropy
%% =========================================================

p = PH(PH > 0);

HH = -sum(p .* log2(p));


%% =========================================================
% Joint entropy
%% =========================================================

p = Pjoint(Pjoint > 0);

HFH = -sum(p .* log2(p));


%% =========================================================
% Mutual information
%
% I(F;H) = H(F) + H(H) - H(F,H)
%% =========================================================

MI = HF + HH - HFH;


%% =========================================================
% Numerical protection
%% =========================================================

MI = max(MI,0);


%% =========================================================
% Normalized Mutual Information
%
% 2I / (H(F) + H(H))
%% =========================================================

if HF + HH > 0

    NMI = 2*MI/(HF+HH);

else

    NMI = NaN;

end


%% =========================================================
% Conditional entropy
%
% H(F|H) = H(F,H) - H(H)
%% =========================================================

H_F_given_H = HFH - HH;


%% =========================================================
% H(H|F)
%% =========================================================

H_H_given_F = HFH - HF;


%% =========================================================
% Numerical protection
%% =========================================================

H_F_given_H = max(H_F_given_H,0);

H_H_given_F = max(H_H_given_F,0);


%% =========================================================
% Normalized conditional entropy
%% =========================================================

if HF > 0

    Normalized_H_F_given_H = ...
        H_F_given_H/HF;

else

    Normalized_H_F_given_H = NaN;

end


if HH > 0

    Normalized_H_H_given_F = ...
        H_H_given_F/HH;

else

    Normalized_H_H_given_F = NaN;

end


%% =========================================================
% Output
%% =========================================================

metrics.HF = HF;

metrics.HH = HH;

metrics.HFH = HFH;

metrics.MI = MI;

metrics.NMI = NMI;

metrics.H_F_given_H = ...
    H_F_given_H;

metrics.H_H_given_F = ...
    H_H_given_F;

metrics.Normalized_H_F_given_H = ...
    Normalized_H_F_given_H;

metrics.Normalized_H_H_given_F = ...
    Normalized_H_H_given_F;

end