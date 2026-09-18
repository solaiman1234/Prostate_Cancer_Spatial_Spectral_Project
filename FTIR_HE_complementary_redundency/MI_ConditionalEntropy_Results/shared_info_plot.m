clc
close all

%% =========================================================
% MULTIMODAL SHARED INFORMATION MAP
%
% X-axis:
% Shared FTIR information with H&E
%
% Y-axis:
% Shared H&E information with FTIR
%
% Color:
% Normalized Mutual Information
%
% Circle   = Non-cancer
% Triangle = Cancer
%% =========================================================


%% =========================================================
% STEP 1
% LOAD MI / CONDITIONAL ENTROPY RESULTS
%% =========================================================

resultFile = ...
'D:\Prostate Cancer Assignment\MI_ConditionalEntropy_Results\data_complementary_redundency.mat';




%% =========================================================
% STEP 2
% CANCER / NON-CANCER FILE LISTS
%% =========================================================

nonCancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


%% =========================================================
% STEP 3
% READ FILENAMES
%% =========================================================

[~,~,nonCancerTable] = xlsread(nonCancerFile);

[~,~,cancerTable] = xlsread(cancerFile);


nonCancerNames = string(nonCancerTable(2:end,1));

cancerNames = string(cancerTable(2:end,1));


%% =========================================================
% STEP 4
% FILENAMES FROM INFORMATION RESULTS
%% =========================================================

resultNames = string(FileName);

% Remove .h5
resultNames = erase(resultNames,'.h5');


%% =========================================================
% STEP 5
% IDENTIFY GROUPS
%% =========================================================

nonCancerIndex = ...
    ismember(resultNames,nonCancerNames);

cancerIndex = ...
    ismember(resultNames,cancerNames);


fprintf('Non-cancer = %d\n',sum(nonCancerIndex));
fprintf('Cancer = %d\n',sum(cancerIndex));


%% =========================================================
% STEP 6
% GET INFORMATION METRICS
%% =========================================================

NMI = NormalizedMI;

FTIR_nonshared = ...
    Normalized_H_FTIR_given_HE;

HE_nonshared = ...
    Normalized_H_HE_given_FTIR;


%% =========================================================
% STEP 7
% CONVERT NON-SHARED -> SHARED
%% =========================================================

FTIR_shared = 1 - FTIR_nonshared;

HE_shared = 1 - HE_nonshared;


%% =========================================================
% STEP 8
% REMOVE INVALID VALUES
%% =========================================================

valid = ...
    isfinite(NMI) & ...
    isfinite(FTIR_shared) & ...
    isfinite(HE_shared);


nonCancerIndex = ...
    nonCancerIndex & valid;

cancerIndex = ...
    cancerIndex & valid;


%% =========================================================
% STEP 9
% SEPARATE GROUPS
%% =========================================================

FTIR_NC = ...
    FTIR_shared(nonCancerIndex);

HE_NC = ...
    HE_shared(nonCancerIndex);

NMI_NC = ...
    NMI(nonCancerIndex);


FTIR_C = ...
    FTIR_shared(cancerIndex);

HE_C = ...
    HE_shared(cancerIndex);

NMI_C = ...
    NMI(cancerIndex);


%% =========================================================
% STEP 10
% GROUP MEDIANS
%% =========================================================

median_FTIR_NC = ...
    median(FTIR_NC,'omitnan');

median_HE_NC = ...
    median(HE_NC,'omitnan');


median_FTIR_C = ...
    median(FTIR_C,'omitnan');

median_HE_C = ...
    median(HE_C,'omitnan');


%% =========================================================
% STEP 11
% CREATE FIGURE
%% =========================================================

figure( ...
    'Color','w', ...
    'Position',[100 100 850 700]);


%% =========================================================
% NON-CANCER
%% =========================================================

scatter( ...
    FTIR_NC, ...
    HE_NC, ...
    70, ...
    NMI_NC, ...
    'o', ...
    'filled', ...
    'MarkerFaceAlpha',0.65, ...
    'MarkerEdgeColor',[0.25 0.25 0.25]);

hold on;


%% =========================================================
% CANCER
%% =========================================================

scatter( ...
    FTIR_C, ...
    HE_C, ...
    85, ...
    NMI_C, ...
    '^', ...
    'filled', ...
    'MarkerFaceAlpha',0.70, ...
    'MarkerEdgeColor',[0.25 0.25 0.25]);


%% =========================================================
% STEP 12
% GROUP MEDIANS
%% =========================================================

plot( ...
    median_FTIR_NC, ...
    median_HE_NC, ...
    'kp', ...
    'MarkerSize',17, ...
    'MarkerFaceColor','w', ...
    'LineWidth',2);


plot( ...
    median_FTIR_C, ...
    median_HE_C, ...
    'kh', ...
    'MarkerSize',17, ...
    'MarkerFaceColor','k', ...
    'LineWidth',2);


%% =========================================================
% STEP 13
% COLORBAR
%% =========================================================

cb = colorbar;

cb.Label.String = ...
    'Normalized Mutual Information';

cb.Label.FontSize = 12;

cb.Label.FontWeight = 'bold';


%% =========================================================
% COLOR LIMITS
%
% Use actual data range instead of 0-1 because your NMI
% values appear to occupy a narrow range.
%% =========================================================

allNMI = [NMI_NC; NMI_C];

cmin = min(allNMI);
cmax = max(allNMI);

if cmax > cmin
    clim([cmin cmax]);
end


%% =========================================================
% STEP 14
% DATA-DRIVEN MEDIAN REFERENCE LINES
%% =========================================================

allFTIR = [FTIR_NC; FTIR_C];
allHE   = [HE_NC; HE_C];

xMedian = median(allFTIR,'omitnan');
yMedian = median(allHE,'omitnan');


xline( ...
    xMedian, ...
    '--', ...
    'LineWidth',1.3);

yline( ...
    yMedian, ...
    '--', ...
    'LineWidth',1.3);


%% =========================================================
% STEP 15
% AXIS LIMITS BASED ON DATA
%% =========================================================

xMin = min(allFTIR);
xMax = max(allFTIR);

yMin = min(allHE);
yMax = max(allHE);


xMargin = 0.05*(xMax-xMin);

yMargin = 0.05*(yMax-yMin);


xlim([ ...
    max(0,xMin-xMargin), ...
    min(1,xMax+xMargin)]);


ylim([ ...
    max(0,yMin-yMargin), ...
    min(1,yMax+yMargin)]);


%% =========================================================
% STEP 16
% LABELS
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
    'Multimodal Shared Information Map', ...
    'FontSize',15, ...
    'FontWeight','bold');


%% =========================================================
% STEP 17
% CLEAN LEGEND
%% =========================================================

h1 = plot( ...
    nan,nan,'o', ...
    'MarkerSize',8, ...
    'MarkerFaceColor',[0.6 0.6 0.6], ...
    'MarkerEdgeColor','k');


h2 = plot( ...
    nan,nan,'^', ...
    'MarkerSize',9, ...
    'MarkerFaceColor',[0.6 0.6 0.6], ...
    'MarkerEdgeColor','k');


h3 = plot( ...
    nan,nan,'kp', ...
    'MarkerSize',12, ...
    'MarkerFaceColor','w', ...
    'LineWidth',1.5);


h4 = plot( ...
    nan,nan,'kh', ...
    'MarkerSize',12, ...
    'MarkerFaceColor','k', ...
    'LineWidth',1.5);


legend( ...
    [h1 h2 h3 h4], ...
    {'Non-cancer', ...
     'Cancer', ...
     'Non-cancer median', ...
     'Cancer median'}, ...
    'Location','best');


%% =========================================================
% STEP 18
% FIGURE FORMATTING
%% =========================================================

grid on;
box on;

set(gca, ...
    'FontSize',12, ...
    'LineWidth',1.3);

axis square;

hold off;


%% =========================================================
% STEP 19
% DISPLAY MEDIAN RESULTS
%% =========================================================

fprintf('\n========================================\n');
fprintf('NON-CANCER\n');
fprintf('========================================\n');

fprintf( ...
    'Median FTIR shared information = %.4f\n', ...
    median_FTIR_NC);

fprintf( ...
    'Median H&E shared information = %.4f\n', ...
    median_HE_NC);

fprintf( ...
    'Median NMI = %.4f\n', ...
    median(NMI_NC,'omitnan'));


fprintf('\n========================================\n');
fprintf('CANCER\n');
fprintf('========================================\n');

fprintf( ...
    'Median FTIR shared information = %.4f\n', ...
    median_FTIR_C);

fprintf( ...
    'Median H&E shared information = %.4f\n', ...
    median_HE_C);

fprintf( ...
    'Median NMI = %.4f\n', ...
    median(NMI_C,'omitnan'));


%% =========================================================
% STEP 20
% CANCER VS NON-CANCER STATISTICS
%% =========================================================

[pFTIR,~,~] = ...
    ranksum(FTIR_NC,FTIR_C);

[pHE,~,~] = ...
    ranksum(HE_NC,HE_C);

[pNMI,~,~] = ...
    ranksum(NMI_NC,NMI_C);


fprintf('\n========================================\n');
fprintf('Cancer vs Non-cancer\n');
fprintf('========================================\n');

fprintf( ...
    'FTIR shared information: p = %.6f\n', ...
    pFTIR);

fprintf( ...
    'H&E shared information: p = %.6f\n', ...
    pHE);

fprintf( ...
    'Normalized MI: p = %.6f\n', ...
    pNMI);

fprintf('\nBonferroni threshold = %.4f\n', ...
    0.05/3);