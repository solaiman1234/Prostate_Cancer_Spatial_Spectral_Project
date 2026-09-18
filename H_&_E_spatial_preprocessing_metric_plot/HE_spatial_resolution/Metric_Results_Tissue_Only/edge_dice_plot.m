clc



%% =========================================================
% H&E CANCER vs NON-CANCER
% FEATURE COMPARISON ACROSS SPATIAL RESOLUTION
%
% Each feature:
% Column 1 = filename
% Column 2 = 1024
% Column 3 = 512
% Column 4 = 256
% Column 5 = 128
% Column 6 = 64
%
% Do NOT use "clear" because feature variables should
% already exist in the MATLAB workspace.
%% =========================================================

clc;
close all;
clearvars -except EdgeDice_data EntropyRetention_data GLCMHomogeneity_data GradientRetention_data HFRetention_data LaplacianRetention_data SSIM_data GLCMContrast_data 

%% =========================================================
% STEP 1
% EXCEL FILES
%% =========================================================

nonCancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


%% =========================================================
% STEP 2
% SPATIAL RESOLUTIONS
%% =========================================================

resolutions = [1024 512 256 128 64];


%% =========================================================
% STEP 3
% DEFINE 8 FEATURES
%% =========================================================

featureData = { ...
    SSIM_data, ...
    EdgeDice_data, ...
    GradientRetention_data, ...
    LaplacianRetention_data, ...
    HFRetention_data, ...
    GLCMContrast_data, ...
    GLCMHomogeneity_data, ...
    EntropyRetention_data};


featureNames = { ...
    'SSIM', ...
    'Edge Dice', ...
    'Gradient Retention', ...
    'Laplacian Retention', ...
    'High-Frequency Retention', ...
    'GLCM Contrast', ...
    'GLCM Homogeneity', ...
    'Entropy Retention'};


%% =========================================================
% STEP 4
% READ CANCER / NON-CANCER EXCEL FILES
%% =========================================================

%cancerTable = readtable( ...
%    cancerFile, ...
%    'VariableNamingRule','preserve');


%nonCancerTable = readtable( ...
%    nonCancerFile, ...
%    'VariableNamingRule','preserve');

[~,~,nonCancerTable] = xlsread(nonCancerFile);
[~,~,cancerTable]    = xlsread(cancerFile);

nonCancerNames = string(nonCancerTable(2:end,1));
cancerNames    = string(cancerTable(2:end,1));

fileNames=EdgeDice_data(:,1);
G_names = erase(fileNames, ".png");
G_names=string(G_names);


%% Find corresponding rows
nonCancerIndex = ismember(G_names, nonCancerNames);
cancerIndex    = ismember(G_names, cancerNames);

fprintf('Number of non-cancer files found = %d\n', ...
    sum(nonCancerIndex));

fprintf('Number of cancer files found = %d\n', ...
    sum(cancerIndex));


%% =========================================================
% STEP 5
% FIND FILENAME COLUMN AUTOMATICALLY
%
% Search column names containing:
% file
% filename
% name
%
% If none are found, use first column.
%% =========================================================

edgedice_noncancer = EdgeDice_data(nonCancerIndex,2:6);
edgedice_noncancer =cell2mat(edgedice_noncancer);
edgedice_cancer    = EdgeDice_data(cancerIndex,2:6);
edgedice_cancer =cell2mat(edgedice_cancer);

clearvars -except edgedice_cancer edgedice_noncancer resolutions

%% Spatial resolutions

resolutions = [1024 512 256 128 64];


fprintf('\n');

fprintf( ...
    'Non-cancer matrix size = %d x %d\n', ...
    size(edgedice_noncancer,1), ...
    size(edgedice_noncancer,2));

fprintf( ...
    'Cancer matrix size = %d x %d\n', ...
    size(edgedice_cancer,1), ...
    size(edgedice_cancer,2));


%% =========================================================
% STEP 8
% CALCULATE MEAN AT EACH RESOLUTION
%% =========================================================

mean_cancer = ...
    mean(edgedice_cancer,1,'omitnan');

mean_noncancer = ...
    mean(edgedice_noncancer,1,'omitnan');


%% =========================================================
% STEP 9
% CALCULATE STANDARD DEVIATION
%% =========================================================

std_cancer = ...
    std(edgedice_cancer,0,1,'omitnan');

std_noncancer = ...
    std(edgedice_noncancer,0,1,'omitnan');


%% =========================================================
% STEP 10
% DISPLAY RESULTS
%% =========================================================

fprintf('\n');
fprintf('Spatial resolutions:\n');

disp([1024 512 256 128 64]);


fprintf('Cancer mean Edge Dice:\n');

disp(mean_cancer);


fprintf('Non-cancer mean Edge Dice:\n');

disp(mean_noncancer);


fprintf('Cancer SD:\n');

disp(std_cancer);


fprintf('Non-cancer SD:\n');

disp(std_noncancer);


%% =========================================================
% STEP 11
% CREATE EQUALLY SPACED X AXIS
%
% DO NOT use:
%
% x = [1024 512 256 128 64]
%
% because MATLAB will use actual numerical spacing.
%% =========================================================

x = 1:5;


resolutionLabels = { ...
    '1024', ...
    '512', ...
    '256', ...
    '128', ...
    '64'};


%% =========================================================
% STEP 12
% PLOT CANCER AND NON-CANCER
%% =========================================================

figure;

errorbar( ...
    x, ...
    mean_cancer, ...
    std_cancer, ...
    '-o', ...
    'LineWidth',2, ...
    'Color','r',...
    'MarkerSize',14, ...
    'CapSize',8);

hold on;


errorbar( ...
    x, ...
    mean_noncancer, ...
    std_noncancer, ...
    '-s', ...
    'LineWidth',2, ...
    'Color','b',...
    'MarkerSize',14, ...
    'CapSize',8);


%% =========================================================
% STEP 13
% X AXIS
%% =========================================================

xticks(x);

xticklabels(resolutionLabels);

xlim([0.7 5.3]);


%% =========================================================
% STEP 14
% Y AXIS
%% =========================================================

ylim([0 1.05]);


%% =========================================================
% STEP 15
% AXIS LABELS
%% =========================================================

xlabel( ...
    'Spatial Resolution', ...
    'FontSize',14, ...
    'FontWeight','bold');


ylabel( ...
    'Edge Dice', ...
    'FontSize',14, ...
    'FontWeight','bold');


%% =========================================================
% STEP 16
% TITLE
%% =========================================================

title( ...
    'Edge Preservation vs Spatial Resolution', ...
    'FontSize',16, ...
    'FontWeight','bold');


%% =========================================================
% STEP 17
% LEGEND
%% =========================================================

legend( ...
    'Cancer', ...
    'Non-Cancer', ...
    'Location','best');


%% =========================================================
% STEP 18
% FIGURE FORMATTING
%% =========================================================

grid on;

box on;


set(gca, ...
    'FontSize',14, ...
    'LineWidth',1.2);


set(gcf, ...
    'Color','w');


hold off;