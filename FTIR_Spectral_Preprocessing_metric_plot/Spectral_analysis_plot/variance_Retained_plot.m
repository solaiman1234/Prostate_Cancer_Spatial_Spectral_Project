clc
close all

%% =========================================================
% G:
% Column 1 = file number
% Column 2 = All / No PCA
% Column 3 = PCA 256
% Column 4 = PCA 64
% Column 5 = PCA 16
% Column 6 = PCA 4
% Column 7 = PCA 1
%% =========================================================

%% Excel files
clearvars -except varianceRetained fileNames

G=[fileNames,varianceRetained];
nonCancerFile = 'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = 'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


[~,~,nonCancerExcel] = xlsread(nonCancerFile);
[~,~,cancerExcel]    = xlsread(cancerFile);

nonCancerNames = string(nonCancerExcel(2:end,1));
cancerNames    = string(cancerExcel(2:end,1));


G_names = erase(G(:,1), ".h5");

nonCancerIndex = ismember(G_names, nonCancerNames);
cancerIndex    = ismember(G_names, cancerNames);

G_nonCancer = G(nonCancerIndex,:);
G_cancer    = G(cancerIndex,:);


fprintf('Number of non-cancer files found = %d\n', ...
    size(G_nonCancer,1));

fprintf('Number of cancer files found = %d\n', ...
    size(G_cancer,1));

variance_nonCancer = str2double(G_nonCancer(:,2:7));
variance_cancer    = str2double(G_cancer(:,2:7));



%% =========================================================
% Calculate group mean
%% =========================================================

mean_nonCancer = mean(variance_nonCancer,1,'omitnan');
mean_cancer    = mean(variance_cancer,1,'omitnan');


%% =========================================================
% Standard deviation
%% =========================================================

std_nonCancer = std(variance_nonCancer,[],1,'omitnan');
std_cancer    = std(variance_cancer,[],1,'omitnan');


%% =========================================================
% 95% confidence interval
%% =========================================================

n_nonCancer = size(variance_nonCancer,1);
n_cancer    = size(variance_cancer,1);

CI_nonCancer = 1.96 * std_nonCancer ./ sqrt(n_nonCancer);
CI_cancer    = 1.96 * std_cancer ./ sqrt(n_cancer);


%% =========================================================
% Display values
%% =========================================================

PCA_names = ["All","256","64","16","4","1"];

ResultTable = table( ...
    PCA_names', ...
    mean_nonCancer', ...
    mean_cancer', ...
    CI_nonCancer', ...
    CI_cancer', ...
    'VariableNames', ...
    {'PCA','NonCancerMean','CancerMean', ...
     'NonCancer95CI','Cancer95CI'});

disp(ResultTable);


%% =========================================================
% Line plot
%% =========================================================

x = 1:6;

figure;

errorbar(x, ...
    mean_nonCancer, ...
    CI_nonCancer, ...
    '-o', ...
    'Color','b',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

hold on;

errorbar(x, ...
    mean_cancer, ...
    CI_cancer, ...
    '-s', ...
    'Color','r',...
    'LineWidth',1.8, ...
    'MarkerSize',14);

xticks(x);
xticklabels(PCA_names);

xlabel('Number of PCA Components','FontSize', 14);
ylabel('Variance Retained (%)','FontSize', 14);

legend('Non-cancer','Cancer', ...
    'Location','best');

title('Variance Retained with Progressive PCA Reduction','FontSize', 18);

grid on;
box on;