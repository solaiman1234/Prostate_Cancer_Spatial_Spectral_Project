clc
close all

%% =========================================================
% REQUIRED VARIABLES ALREADY IN WORKSPACE:
%
% meanProcessedSpectrum = 227 x 6 x 965
% fileNames             = 227 x 1 string
%
% PCA order:
% 1 = All
% 2 = 256
% 3 = 64
% 4 = 16
% 5 = 4
% 6 = 1
%% =========================================================

clearvars -except meanProcessedSpectrum fileNames


%% =========================================================
% Excel files containing cancer / non-cancer core names
%% =========================================================

nonCancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


%% =========================================================
% Read Excel files
%% =========================================================

[~,~,nonCancerExcel] = xlsread(nonCancerFile);
[~,~,cancerExcel]    = xlsread(cancerFile);


%% =========================================================
% Get core names from first column
% First row assumed to be header
%% =========================================================

nonCancerNames = string(nonCancerExcel(2:end,1));
cancerNames    = string(cancerExcel(2:end,1));


%% Remove empty entries
nonCancerNames = nonCancerNames(nonCancerNames ~= "");
cancerNames    = cancerNames(cancerNames ~= "");


%% Standardize names
nonCancerNames = lower(strtrim(nonCancerNames));
cancerNames    = lower(strtrim(cancerNames));


%% =========================================================
% File names from MATLAB results
%
% fileNames contains .h5
% Excel names do not contain .h5
%% =========================================================

G_names = lower(strtrim(erase(fileNames, ".h5")));


%% =========================================================
% Find cancer and non-cancer rows
%% =========================================================

nonCancerIndex = ismember(G_names, nonCancerNames);
cancerIndex    = ismember(G_names, cancerNames);


%% Check numbers
fprintf('Number of non-cancer files = %d\n', ...
    sum(nonCancerIndex));

fprintf('Number of cancer files     = %d\n', ...
    sum(cancerIndex));


%% =========================================================
% PCA labels
%% =========================================================

PCA_names = ["All","256","64","16","4","1"];


%% =========================================================
% Original spectral information
%
% Original = 1479 channels
%
% Removed noisy channels:
%% =========================================================

idx_remove = [ ...
    1:12, ...
    204:281, ...
    676:779, ...
    909:1064, ...
    1316:1479];


%% Channels that were retained
idx_keep = setdiff(1:1479, idx_remove);


fprintf('Number of retained channels = %d\n', ...
    length(idx_keep));

% Should print:
% 965


%% =========================================================
% Define continuous retained spectral regions
%
% These are ORIGINAL channel numbers
%% =========================================================

segments = [ ...
      13   203;
     282   675;
     780   908;
    1065  1315];


%% =========================================================
% Create figure
%% =========================================================

figure('Position',[100 100 1600 850]);

t = tiledlayout(2,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% =========================================================
% LOOP THROUGH 6 PCA LEVELS
%% =========================================================

for k = 1:6

    %% -----------------------------------------------------
    % Extract cancer spectra
    %
    % Expected:
    % 127 x 965
    %% -----------------------------------------------------

    cancerData = squeeze( ...
        meanProcessedSpectrum(cancerIndex,k,:));


    %% -----------------------------------------------------
    % Extract non-cancer spectra
    %
    % Expected:
    % 100 x 965
    %% -----------------------------------------------------

    nonCancerData = squeeze( ...
        meanProcessedSpectrum(nonCancerIndex,k,:));


    %% =====================================================
    % Mean spectrum
    %% =====================================================

    meanCancer = mean( ...
        cancerData,1,'omitnan');

    meanNonCancer = mean( ...
        nonCancerData,1,'omitnan');


    %% =====================================================
    % Standard deviation
    %% =====================================================

    stdCancer = std( ...
        cancerData,[],1,'omitnan');

    stdNonCancer = std( ...
        nonCancerData,[],1,'omitnan');


    %% =====================================================
    % Number of valid observations at each channel
    %% =====================================================

    nCancerChannel = ...
        sum(~isnan(cancerData),1);

    nNonCancerChannel = ...
        sum(~isnan(nonCancerData),1);


    %% =====================================================
    % Standard error
    %% =====================================================

    seCancer = ...
        stdCancer ./ sqrt(nCancerChannel);

    seNonCancer = ...
        stdNonCancer ./ sqrt(nNonCancerChannel);


    %% =====================================================
    % Approximate 95% confidence intervals
    %% =====================================================

    CI_Cancer = ...
        1.96 * seCancer;

    CI_NonCancer = ...
        1.96 * seNonCancer;


    %% =====================================================
    % Restore 965 channels into original 1479 locations
    %
    % Removed channels remain NaN
    %% =====================================================

    meanCancer_full = nan(1,1479);

    meanNonCancer_full = nan(1,1479);

    CI_Cancer_full = nan(1,1479);

    CI_NonCancer_full = nan(1,1479);


    %% Insert retained values
    meanCancer_full(idx_keep) = ...
        meanCancer;

    meanNonCancer_full(idx_keep) = ...
        meanNonCancer;

    CI_Cancer_full(idx_keep) = ...
        CI_Cancer;

    CI_NonCancer_full(idx_keep) = ...
        CI_NonCancer;


    %% =====================================================
    % Plot
    %% =====================================================

    nexttile;

    hold on;


    %% =====================================================
    % Plot 95% confidence intervals separately for
    % each continuous spectral region
    %% =====================================================

    for s = 1:size(segments,1)

        xx = ...
            segments(s,1):segments(s,2);


        %% -------------------------------------------------
        % Non-cancer confidence interval
        %% -------------------------------------------------

        hNCci = fill( ...
            [xx fliplr(xx)], ...
            [meanNonCancer_full(xx) + ...
             CI_NonCancer_full(xx), ...
             fliplr( ...
             meanNonCancer_full(xx) - ...
             CI_NonCancer_full(xx))], ...
            [0.80 0.80 0.80], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.30);


        %% -------------------------------------------------
        % Cancer confidence interval
        %% -------------------------------------------------

        hCci = fill( ...
            [xx fliplr(xx)], ...
            [meanCancer_full(xx) + ...
             CI_Cancer_full(xx), ...
             fliplr( ...
             meanCancer_full(xx) - ...
             CI_Cancer_full(xx))], ...
            [0.65 0.65 0.65], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.30);


        %% Do not repeat CI entries in legend
        if s > 1

            hNCci.HandleVisibility = 'off';
            hCci.HandleVisibility  = 'off';

        end

    end


    %% =====================================================
    % Plot mean spectra
    %
    % NaN regions automatically appear as gaps
    %% =====================================================

    hNC = plot( ...
        1:1479, ...
        meanNonCancer_full, ...
         'Color','b',...
        'LineWidth',1.5);


    hC = plot( ...
        1:1479, ...
        meanCancer_full, ...
         'Color','r',...
        'LineWidth',1.5);


    %% =====================================================
    % Axis formatting
    %% =====================================================

    xlim([1 1479]);

    xlabel('Original Spectral Channel','FontSize',14);

    ylabel('2nd Derivative','FontSize', 14);

    title(['PCA: ',char(PCA_names(k))],'FontSize',16);

    grid on;
    box on;


    %% =====================================================
    % Legend only in first panel
    %% =====================================================

    if k == 1

        legend( ...
            [hNCci hCci hNC hC], ...
            {'Non-cancer 95% CI', ...
             'Cancer 95% CI', ...
             'Non-cancer mean', ...
             'Cancer mean'}, ...
            'Location','best');

    end

end


%% =========================================================
% Overall title
%% =========================================================

title(t, ...
    'Cancer vs Non-cancer Mean Processed FTIR Spectra','FontSize', 16);