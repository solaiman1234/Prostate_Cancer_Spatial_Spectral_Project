clc
close all

%% Keep required variables
clearvars -except NRMSE meanSAM fileNames

G=[fileNames,NRMSE,meanSAM];

%% =========================================================
% Excel files
%% =========================================================

nonCancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


[~,~,nonCancerExcel] = xlsread(nonCancerFile);
[~,~,cancerExcel]    = xlsread(cancerFile);


%% =========================================================
% Get cancer / non-cancer file names
%% =========================================================

nonCancerNames = string(nonCancerExcel(2:end,1));
cancerNames    = string(cancerExcel(2:end,1));


%% Remove .h5 from combined file names
G_names = erase(fileNames, ".h5");


%% Find corresponding rows
nonCancerIndex = ismember(G_names, nonCancerNames);
cancerIndex    = ismember(G_names, cancerNames);


fprintf('Number of non-cancer files found = %d\n', ...
    sum(nonCancerIndex));

fprintf('Number of cancer files found = %d\n', ...
    sum(cancerIndex));


%% =========================================================
% Separate NRMSE
%% =========================================================

NRMSE_nonCancer = NRMSE(nonCancerIndex,:);
NRMSE_cancer    = NRMSE(cancerIndex,:);


%% =========================================================
% Separate Mean SAM
%% =========================================================

SAM_nonCancer = meanSAM(nonCancerIndex,:);
SAM_cancer    = meanSAM(cancerIndex,:);


%% =========================================================
% Calculate means
%% =========================================================

mean_NRMSE_nonCancer = mean(NRMSE_nonCancer,1,'omitnan');
mean_NRMSE_cancer    = mean(NRMSE_cancer,1,'omitnan');

mean_SAM_nonCancer = mean(SAM_nonCancer,1,'omitnan');
mean_SAM_cancer    = mean(SAM_cancer,1,'omitnan');


%% =========================================================
% Standard deviations
%% =========================================================

std_NRMSE_nonCancer = std(NRMSE_nonCancer,[],1,'omitnan');
std_NRMSE_cancer    = std(NRMSE_cancer,[],1,'omitnan');

std_SAM_nonCancer = std(SAM_nonCancer,[],1,'omitnan');
std_SAM_cancer    = std(SAM_cancer,[],1,'omitnan');


%% =========================================================
% 95% confidence intervals
%% =========================================================

n_nonCancer = sum(nonCancerIndex);
n_cancer    = sum(cancerIndex);


CI_NRMSE_nonCancer = ...
    1.96 * std_NRMSE_nonCancer ./ sqrt(n_nonCancer);

CI_NRMSE_cancer = ...
    1.96 * std_NRMSE_cancer ./ sqrt(n_cancer);


CI_SAM_nonCancer = ...
    1.96 * std_SAM_nonCancer ./ sqrt(n_nonCancer);

CI_SAM_cancer = ...
    1.96 * std_SAM_cancer ./ sqrt(n_cancer);


%% =========================================================
% PCA labels
%% =========================================================

PCA_names = ["All","256","64","16","4","1"];

x = 1:6;


%% =========================================================
% FIGURE 2
% A = NRMSE
% B = Mean SAM
%% =========================================================

figure;

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');


%% ---------------------------------------------------------
% A. NRMSE
%% ---------------------------------------------------------

nexttile;

errorbar(x, ...
    mean_NRMSE_nonCancer, ...
    CI_NRMSE_nonCancer, ...
    '-o', ...
    'Color','b',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

hold on;

errorbar(x, ...
    mean_NRMSE_cancer, ...
    CI_NRMSE_cancer, ...
    '-s', ...
    'Color','r',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

xticks(x);
xticklabels(PCA_names);

xlabel('Number of PCA Components','FontSize', 14);
ylabel('NRMSE','FontSize', 14);

title('(A) Reconstruction Error','FontSize', 18);

legend('Non-cancer','Cancer','Location','best');

grid on;
box on;


%% ---------------------------------------------------------
% B. Mean SAM
%% ---------------------------------------------------------

nexttile;

errorbar(x, ...
    mean_SAM_nonCancer, ...
    CI_SAM_nonCancer, ...
    '-o', ...
     'Color','b',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

hold on;

errorbar(x, ...
    mean_SAM_cancer, ...
    CI_SAM_cancer, ...
    '-s', ...
     'Color','r',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

xticks(x);
xticklabels(PCA_names);

xlabel('Number of PCA Components','FontSize', 14);
ylabel('Mean SAM (degrees)','FontSize', 14);

title('(B) Spectral Angle','FontSize', 16);

legend('Non-cancer','Cancer','Location','best');

grid on;
box on;


sgtitle('Spectral Reconstruction Quality with Progressive PCA Reduction','FontSize', 16);