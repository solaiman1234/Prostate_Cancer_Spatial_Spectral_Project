clc
close all

%% Keep required variables
clearvars -except derivativeCorrelation fileNames

G=[fileNames,derivativeCorrelation];
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
% Separate derivative correlation
%% =========================================================

DC_nonCancer = derivativeCorrelation(nonCancerIndex,:);
DC_cancer    = derivativeCorrelation(cancerIndex,:);


%% =========================================================
% Calculate group mean
%% =========================================================

mean_nonCancer = mean(DC_nonCancer,1,'omitnan');
mean_cancer    = mean(DC_cancer,1,'omitnan');


%% =========================================================
% Standard deviation
%% =========================================================

std_nonCancer = std(DC_nonCancer,[],1,'omitnan');
std_cancer    = std(DC_cancer,[],1,'omitnan');


%% =========================================================
% 95% confidence interval
%% =========================================================

n_nonCancer = size(DC_nonCancer,1);
n_cancer    = size(DC_cancer,1);

CI_nonCancer = ...
    1.96 * std_nonCancer ./ sqrt(n_nonCancer);

CI_cancer = ...
    1.96 * std_cancer ./ sqrt(n_cancer);


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

ylabel('Derivative Correlation','FontSize', 14);

legend('Non-cancer','Cancer', ...
    'Location','best');

title('Preservation of SG Second-Derivative Spectra','FontSize', 16);

grid on;
box on;