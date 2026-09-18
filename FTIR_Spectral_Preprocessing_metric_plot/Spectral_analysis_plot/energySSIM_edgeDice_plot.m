clc
close all

%% Keep required variables
clearvars -except energySSIM edgeDice fileNames

G=[fileNames,energySSIM,edgeDice];
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


%% Remove .h5 from combined filenames
G_names = erase(fileNames, ".h5");


%% Find corresponding rows
nonCancerIndex = ismember(G_names, nonCancerNames);
cancerIndex    = ismember(G_names, cancerNames);


fprintf('Number of non-cancer files found = %d\n', ...
    sum(nonCancerIndex));

fprintf('Number of cancer files found = %d\n', ...
    sum(cancerIndex));


%% =========================================================
% Separate Energy-map SSIM
%% =========================================================

SSIM_nonCancer = energySSIM(nonCancerIndex,:);
SSIM_cancer    = energySSIM(cancerIndex,:);


%% =========================================================
% Separate Edge Dice
%% =========================================================

Dice_nonCancer = edgeDice(nonCancerIndex,:);
Dice_cancer    = edgeDice(cancerIndex,:);


%% =========================================================
% Calculate means
%% =========================================================

mean_SSIM_nonCancer = mean(SSIM_nonCancer,1,'omitnan');
mean_SSIM_cancer    = mean(SSIM_cancer,1,'omitnan');

mean_Dice_nonCancer = mean(Dice_nonCancer,1,'omitnan');
mean_Dice_cancer    = mean(Dice_cancer,1,'omitnan');


%% =========================================================
% Standard deviations
%% =========================================================

std_SSIM_nonCancer = std(SSIM_nonCancer,[],1,'omitnan');
std_SSIM_cancer    = std(SSIM_cancer,[],1,'omitnan');

std_Dice_nonCancer = std(Dice_nonCancer,[],1,'omitnan');
std_Dice_cancer    = std(Dice_cancer,[],1,'omitnan');


%% =========================================================
% 95% confidence intervals
%% =========================================================

n_nonCancer = sum(nonCancerIndex);
n_cancer    = sum(cancerIndex);


CI_SSIM_nonCancer = ...
    1.96 * std_SSIM_nonCancer ./ sqrt(n_nonCancer);

CI_SSIM_cancer = ...
    1.96 * std_SSIM_cancer ./ sqrt(n_cancer);


CI_Dice_nonCancer = ...
    1.96 * std_Dice_nonCancer ./ sqrt(n_nonCancer);

CI_Dice_cancer = ...
    1.96 * std_Dice_cancer ./ sqrt(n_cancer);


%% =========================================================
% PCA labels
%% =========================================================

PCA_names = ["All","256","64","16","4","1"];

x = 1:6;


%% =========================================================
% FIGURE 4
% A = Energy-map SSIM
% B = Edge Dice
%% =========================================================

figure;

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');


%% ---------------------------------------------------------
% A. Energy-map SSIM
%% ---------------------------------------------------------

nexttile;

errorbar(x, ...
    mean_SSIM_nonCancer, ...
    CI_SSIM_nonCancer, ...
    '-o', ...
     'Color','b',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

hold on;

errorbar(x, ...
    mean_SSIM_cancer, ...
    CI_SSIM_cancer, ...
    '-s', ...
     'Color','r',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

xticks(x);
xticklabels(PCA_names);

xlabel('Number of PCA Components','FontSize', 14);
ylabel('SSIM','FontSize', 14);

title('(A) Spectral Energy-map SSIM','FontSize', 16);

legend('Non-cancer','Cancer', ...
    'Location','best');

ylim([0 1]);

grid on;
box on;


%% ---------------------------------------------------------
% B. Edge Dice
%% ---------------------------------------------------------

nexttile;

errorbar(x, ...
    mean_Dice_nonCancer, ...
    CI_Dice_nonCancer, ...
    '-o', ...
     'Color','b',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

hold on;

errorbar(x, ...
    mean_Dice_cancer, ...
    CI_Dice_cancer, ...
    '-s', ...
     'Color','r',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

xticks(x);
xticklabels(PCA_names);

xlabel('Number of PCA Components','FontSize', 14);
ylabel('Edge Dice','FontSize', 14);

title('(B) Spatial Boundary Preservation','FontSize', 16);

legend('Non-cancer','Cancer', ...
    'Location','best');

ylim([0 1]);

grid on;
box on;


sgtitle('Spatial Preservation with Progressive PCA Reduction','FontSize', 16);