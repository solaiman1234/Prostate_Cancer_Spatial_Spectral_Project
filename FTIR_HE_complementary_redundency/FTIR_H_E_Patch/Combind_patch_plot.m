clc
%close all
clear

%% =========================================================
% LOAD PATCH-LEVEL MULTIVARIATE RESULTS
%% =========================================================

load( ...
'D:\Prostate Cancer Assignment\Patch_Multivariate_Information\Patch_Multivariate_Information.mat');


%% =========================================================
% REMOVE INVALID VALUES
%% =========================================================

valid = ...
    isfinite(Results.NormalizedMI) & ...
    isfinite(Results.SharedFTIR) & ...
    isfinite(Results.SharedHE);

T = Results(valid,:);


%% =========================================================
% DEFINE GROUPS
%% =========================================================

NC = T.Group == "Non-cancer";
C  = T.Group == "Cancer";


fprintf('Number of cancer samples     = %d\n',sum(C));
fprintf('Number of non-cancer samples = %d\n',sum(NC));


%% =========================================================
% CALCULATE GROUP MEDIANS
%% =========================================================

medianNMI_NC = median(T.NormalizedMI(NC),'omitnan');
medianNMI_C  = median(T.NormalizedMI(C),'omitnan');

medianFTIR_NC = median(T.SharedFTIR(NC),'omitnan');
medianFTIR_C  = median(T.SharedFTIR(C),'omitnan');

medianHE_NC = median(T.SharedHE(NC),'omitnan');
medianHE_C  = median(T.SharedHE(C),'omitnan');


%% =========================================================
% DISPLAY MEDIANS
%% =========================================================

fprintf('\n============================================\n');
fprintf('NON-CANCER MEDIAN VALUES\n');
fprintf('============================================\n');

fprintf('NMI         = %.4f\n',medianNMI_NC);
fprintf('Shared FTIR = %.4f\n',medianFTIR_NC);
fprintf('Shared H&E  = %.4f\n',medianHE_NC);


fprintf('\n============================================\n');
fprintf('CANCER MEDIAN VALUES\n');
fprintf('============================================\n');

fprintf('NMI         = %.4f\n',medianNMI_C);
fprintf('Shared FTIR = %.4f\n',medianFTIR_C);
fprintf('Shared H&E  = %.4f\n',medianHE_C);



%% =========================================================
% FIGURE 1
%
% PATCH-LEVEL MULTIVARIATE SHARED INFORMATION MAP
%% =========================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 900 750]);


%% ---------------------------------------------------------
% NON-CANCER
%
% Circle = non-cancer
%% ---------------------------------------------------------

hNC = scatter( ...
    T.SharedFTIR(NC), ...
    T.SharedHE(NC), ...
    65, ...
    T.NormalizedMI(NC), ...
    'o', ...
    'filled', ...
    'MarkerFaceAlpha',0.70, ...
    'MarkerEdgeColor',[0.2 0.2 0.2]);

hold on;


%% ---------------------------------------------------------
% CANCER
%
% Triangle = cancer
%% ---------------------------------------------------------

hC = scatter( ...
    T.SharedFTIR(C), ...
    T.SharedHE(C), ...
    80, ...
    T.NormalizedMI(C), ...
    '^', ...
    'filled', ...
    'MarkerFaceAlpha',0.75, ...
    'MarkerEdgeColor',[0.2 0.2 0.2]);


%% =========================================================
% NON-CANCER MEDIAN
%% =========================================================

hNCmed = plot( ...
    medianFTIR_NC, ...
    medianHE_NC, ...
    'kp', ...
    'MarkerSize',17, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',2);


%% =========================================================
% CANCER MEDIAN
%% =========================================================

hCmed = plot( ...
    medianFTIR_C, ...
    medianHE_C, ...
    'kh', ...
    'MarkerSize',17, ...
    'MarkerFaceColor','k', ...
    'MarkerEdgeColor','k', ...
    'LineWidth',2);


%% =========================================================
% GLOBAL MEDIAN REFERENCE LINES
%
% These are descriptive reference lines only.
%% =========================================================

overallFTIRmedian = ...
    median(T.SharedFTIR,'omitnan');

overallHEmedian = ...
    median(T.SharedHE,'omitnan');


xline( ...
    overallFTIRmedian, ...
    '--', ...
    'LineWidth',1.4, ...
    'Color',[0.4 0.4 0.4], ...
    'HandleVisibility','off');


yline( ...
    overallHEmedian, ...
    '--', ...
    'LineWidth',1.4, ...
    'Color',[0.4 0.4 0.4], ...
    'HandleVisibility','off');


%% =========================================================
% COLORBAR
%
% Color = Normalized Mutual Information
%% =========================================================

cb = colorbar;

cb.Label.String = ...
    'Normalized Mutual Information';

cb.Label.FontSize = 12;

cb.Label.FontWeight = 'bold';


%% =========================================================
% COLOR RANGE
%% =========================================================

NMImin = min(T.NormalizedMI);
NMImax = max(T.NormalizedMI);

if NMImax > NMImin

    clim([NMImin NMImax]);

end

colormap(parula);


%% =========================================================
% LEGEND
%% =========================================================

legend( ...
    [hNC hC hNCmed hCmed], ...
    { ...
    'Non-cancer', ...
    'Cancer', ...
    'Non-cancer median', ...
    'Cancer median'}, ...
    'Location','best', ...
    'FontSize',10);


%% =========================================================
% AXIS LABELS
%% =========================================================

xlabel( ...
    'Shared FTIR Information with H&E', ...
    'FontSize',13, ...
    'FontWeight','bold');


ylabel( ...
    'Shared H&E Information with FTIR', ...
    'FontSize',13, ...
    'FontWeight','bold');


title( ...
    'Patch-level Multivariate Shared Information Map', ...
    'FontSize',14, ...
    'FontWeight','bold');


%% =========================================================
% FORMATTING
%% =========================================================

grid on
box on

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.2);

axis square;



%% =========================================================
% FIGURE 2
%
% INDIVIDUAL DISTRIBUTIONS
%% =========================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 1250 430]);


%% =========================================================
% FORCE CATEGORY ORDER
%
% Cancer first
% Non-cancer second
%% =========================================================

groupCat = categorical( ...
    T.Group, ...
    ["Cancer","Non-cancer"], ...
    ["Cancer","Non-cancer"]);



%% =========================================================
% PANEL 1
% NORMALIZED MUTUAL INFORMATION
%% =========================================================

subplot(1,3,1)


boxchart( ...
    groupCat, ...
    T.NormalizedMI);

hold on


swarmchart( ...
    groupCat, ...
    T.NormalizedMI, ...
    18, ...
    'filled', ...
    'MarkerFaceAlpha',0.75);


ylabel( ...
    'Normalized Mutual Information', ...
    'FontWeight','bold');


xlabel('');


title( ...
    'Cross-modal NMI', ...
    'FontWeight','bold');


grid on
box on


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.1);



%% =========================================================
% PANEL 2
% FTIR INFORMATION SHARED WITH H&E
%% =========================================================

subplot(1,3,2)


boxchart( ...
    groupCat, ...
    T.SharedFTIR);

hold on


swarmchart( ...
    groupCat, ...
    T.SharedFTIR, ...
    18, ...
    'filled', ...
    'MarkerFaceAlpha',0.75);


ylabel( ...
    'Shared Information Fraction', ...
    'FontWeight','bold');


xlabel('');


title( ...
    'FTIR Shared with H&E', ...
    'FontWeight','bold');


grid on
box on


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.1);



%% =========================================================
% PANEL 3
% H&E INFORMATION SHARED WITH FTIR
%% =========================================================

subplot(1,3,3)


boxchart( ...
    groupCat, ...
    T.SharedHE);

hold on


swarmchart( ...
    groupCat, ...
    T.SharedHE, ...
    18, ...
    'filled', ...
    'MarkerFaceAlpha',0.75);


ylabel( ...
    'Shared Information Fraction', ...
    'FontWeight','bold');


xlabel('');


title( ...
    'H&E Shared with FTIR', ...
    'FontWeight','bold');


grid on
box on


set(gca, ...
    'FontSize',11, ...
    'LineWidth',1.1);



%% =========================================================
% STATISTICAL COMPARISON
%
% Wilcoxon rank-sum test:
% Cancer vs Non-cancer
%% =========================================================

pNMI = ranksum( ...
    T.NormalizedMI(C), ...
    T.NormalizedMI(NC));


pFTIR = ranksum( ...
    T.SharedFTIR(C), ...
    T.SharedFTIR(NC));


pHE = ranksum( ...
    T.SharedHE(C), ...
    T.SharedHE(NC));


%% =========================================================
% DISPLAY STATISTICS
%% =========================================================

fprintf('\n============================================\n');
fprintf('CANCER VS NON-CANCER STATISTICS\n');
fprintf('============================================\n');

fprintf( ...
    'Normalized MI:       p = %.6f\n', ...
    pNMI);

fprintf( ...
    'FTIR shared with H&E: p = %.6f\n', ...
    pFTIR);

fprintf( ...
    'H&E shared with FTIR: p = %.6f\n', ...
    pHE);


%% =========================================================
% BONFERRONI CORRECTION
%% =========================================================

alphaCorrected = 0.05/3;

fprintf('\nBonferroni corrected alpha = %.4f\n', ...
    alphaCorrected);


%% =========================================================
% SIGNIFICANCE MESSAGE
%% =========================================================

fprintf('\n');

if pNMI < alphaCorrected
    fprintf('NMI difference is statistically significant.\n');
else
    fprintf('NMI difference is not statistically significant.\n');
end


if pFTIR < alphaCorrected
    fprintf('FTIR shared-information difference is statistically significant.\n');
else
    fprintf('FTIR shared-information difference is not statistically significant.\n');
end


if pHE < alphaCorrected
    fprintf('H&E shared-information difference is statistically significant.\n');
else
    fprintf('H&E shared-information difference is not statistically significant.\n');
end