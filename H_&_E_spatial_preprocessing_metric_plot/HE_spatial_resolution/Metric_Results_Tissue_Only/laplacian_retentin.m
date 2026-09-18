clc;
close all;

%% =========================================================
% H&E CANCER vs NON-CANCER
% LAPLACIAN RETENTION COMPARISON ACROSS SPATIAL RESOLUTION
%
% LaplacianRetention_data:
% Column 1 = filename
% Column 2 = 1024
% Column 3 = 512
% Column 4 = 256
% Column 5 = 128
% Column 6 = 64
%% =========================================================


%% =========================================================
% KEEP REQUIRED VARIABLES
%% =========================================================

clearvars -except EdgeDice_data EntropyRetention_data ...
    GLCMHomogeneity_data GradientRetention_data ...
    HFRetention_data LaplacianRetention_data ...
    SSIM_data GLCMContrast_data


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
% READ CANCER / NON-CANCER EXCEL FILES
%% =========================================================

[~,~,nonCancerTable] = xlsread(nonCancerFile);

[~,~,cancerTable] = xlsread(cancerFile);


nonCancerNames = ...
    string(nonCancerTable(2:end,1));

cancerNames = ...
    string(cancerTable(2:end,1));


%% =========================================================
% STEP 4
% GET FILENAMES FROM LAPLACIAN RETENTION DATA
%% =========================================================

fileNames = ...
    LaplacianRetention_data(:,1);


G_names = ...
    erase(fileNames,'.png');


G_names = ...
    string(G_names);


%% =========================================================
% STEP 5
% FIND CORRESPONDING ROWS
%% =========================================================

nonCancerIndex = ...
    ismember(G_names,nonCancerNames);

cancerIndex = ...
    ismember(G_names,cancerNames);


fprintf( ...
    'Number of non-cancer files found = %d\n', ...
    sum(nonCancerIndex));


fprintf( ...
    'Number of cancer files found = %d\n', ...
    sum(cancerIndex));


%% =========================================================
% STEP 6
% SEPARATE NON-CANCER AND CANCER
%% =========================================================

laplacianretention_noncancer = ...
    LaplacianRetention_data(nonCancerIndex,2:6);


laplacianretention_noncancer = ...
    cell2mat(laplacianretention_noncancer);


laplacianretention_cancer = ...
    LaplacianRetention_data(cancerIndex,2:6);


laplacianretention_cancer = ...
    cell2mat(laplacianretention_cancer);


%% =========================================================
% KEEP ONLY REQUIRED VARIABLES
%% =========================================================

clearvars -except laplacianretention_cancer ...
    laplacianretention_noncancer resolutions


%% =========================================================
% SPATIAL RESOLUTIONS
%% =========================================================

resolutions = [1024 512 256 128 64];


%% =========================================================
% CHECK MATRIX SIZE
%% =========================================================

fprintf('\n');


fprintf( ...
    'Non-cancer matrix size = %d x %d\n', ...
    size(laplacianretention_noncancer,1), ...
    size(laplacianretention_noncancer,2));


fprintf( ...
    'Cancer matrix size = %d x %d\n', ...
    size(laplacianretention_cancer,1), ...
    size(laplacianretention_cancer,2));


%% =========================================================
% STEP 7
% CALCULATE MEAN AT EACH RESOLUTION
%% =========================================================

mean_cancer = ...
    mean(laplacianretention_cancer,1,'omitnan');


mean_noncancer = ...
    mean(laplacianretention_noncancer,1,'omitnan');


%% =========================================================
% STEP 8
% CALCULATE STANDARD DEVIATION
%% =========================================================

std_cancer = ...
    std(laplacianretention_cancer,0,1,'omitnan');


std_noncancer = ...
    std(laplacianretention_noncancer,0,1,'omitnan');


%% =========================================================
% STEP 9
% DISPLAY RESULTS
%% =========================================================

fprintf('\n');

fprintf('Spatial resolutions:\n');

disp([1024 512 256 128 64]);


fprintf('Cancer mean Laplacian Retention:\n');

disp(mean_cancer);


fprintf('Non-cancer mean Laplacian Retention:\n');

disp(mean_noncancer);


fprintf('Cancer SD:\n');

disp(std_cancer);


fprintf('Non-cancer SD:\n');

disp(std_noncancer);


%% =========================================================
% STEP 10
% EQUALLY SPACED X AXIS
%% =========================================================

x = 1:5;


resolutionLabels = { ...
    '1024', ...
    '512', ...
    '256', ...
    '128', ...
    '64'};


%% =========================================================
% STEP 11
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
% STEP 12
% X AXIS
%% =========================================================

xticks(x);

xticklabels(resolutionLabels);

xlim([0.7 5.3]);


%% =========================================================
% STEP 13
% Y AXIS
%% =========================================================

ylim([0 1.05]);


%% =========================================================
% STEP 14
% AXIS LABELS
%% =========================================================

xlabel( ...
    'Spatial Resolution', ...
    'FontSize',14, ...
    'FontWeight','bold');


ylabel( ...
    'Laplacian Retention', ...
    'FontSize',14, ...
    'FontWeight','bold');


%% =========================================================
% STEP 15
% TITLE
%% =========================================================

title( ...
    'Laplacian Retention vs Spatial Resolution', ...
    'FontSize',14, ...
    'FontWeight','bold');


%% =========================================================
% STEP 16
% LEGEND
%% =========================================================

legend( ...
    'Cancer', ...
    'Non-Cancer', ...
    'Location','best');


%% =========================================================
% STEP 17
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