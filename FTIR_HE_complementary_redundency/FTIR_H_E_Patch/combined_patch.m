clc
clear
close all

%% =========================================================
% INPUT FOLDERS
%% =========================================================

ftirFolder = ...
'D:\Prostate Cancer Assignment\spectral_annotation_965';

heFolder = ...
'D:\Prostate Cancer Assignment\aligned_he_annotation_only';


nonCancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_non_cancer.xlsx';

cancerFile = ...
'D:\Prostate Cancer Assignment\master_sheet_annotation_only_cancer.xlsx';


outputFolder = ...
'D:\Prostate Cancer Assignment\Patch_Multivariate_Information';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end


%% =========================================================
% PARAMETERS
%% =========================================================

patchSize = 16;

FTIRsmoothPC = 64;

patchFeaturePC = 5;

minimumTissueFraction = 0.80;

datasetName = '/spectra';


%% =========================================================
% READ CANCER / NON-CANCER NAMES
%% =========================================================
%% =========================================================
% READ CANCER / NON-CANCER FILE NAMES
%% =========================================================

[~,~,nonCancerTable] = xlsread(nonCancerFile);
[~,~,cancerTable]    = xlsread(cancerFile);

% First column contains image/file names
% Skip first row because it is the header
nonCancerNames = string(nonCancerTable(2:end,1));
cancerNames    = string(cancerTable(2:end,1));


%% =========================================================
% CLEAN FILENAMES
%% =========================================================

nonCancerNames = strip(nonCancerNames);
cancerNames    = strip(cancerNames);

nonCancerNames = lower(nonCancerNames);
cancerNames    = lower(cancerNames);

% Remove possible extensions
nonCancerNames = erase(nonCancerNames,'.png');
nonCancerNames = erase(nonCancerNames,'.h5');

cancerNames = erase(cancerNames,'.png');
cancerNames = erase(cancerNames,'.h5');


%% =========================================================
% REMOVE EMPTY ENTRIES
%% =========================================================

nonCancerNames = ...
    nonCancerNames( ...
    ~ismissing(nonCancerNames) & ...
    strlength(nonCancerNames) > 0);

cancerNames = ...
    cancerNames( ...
    ~ismissing(cancerNames) & ...
    strlength(cancerNames) > 0);


%% =========================================================
% CHECK
%% =========================================================

fprintf('Non-cancer names loaded = %d\n',length(nonCancerNames));
fprintf('Cancer names loaded     = %d\n',length(cancerNames));

nonCancerNames = lower(strip(nonCancerNames));
cancerNames    = lower(strip(cancerNames));

nonCancerNames = erase(nonCancerNames,{'.png','.h5'});
cancerNames    = erase(cancerNames,{'.png','.h5'});

ftirFiles = dir(fullfile(ftirFolder,'*.h5'));

nFiles = length(ftirFiles);

fprintf('Total FTIR files = %d\n',nFiles);


FileName = strings(nFiles,1);

Group = strings(nFiles,1);

NumberPatches = nan(nFiles,1);

MI = nan(nFiles,1);

EntropyFTIR = nan(nFiles,1);

EntropyHE = nan(nFiles,1);

NormalizedMI = nan(nFiles,1);

SharedFTIR = nan(nFiles,1);

SharedHE = nan(nFiles,1);

%% =========================================================
% PROCESS ALL IMAGE PAIRS
%% =========================================================

for i = 1:nFiles

    fprintf('\n========================================\n');
    fprintf('Processing %d / %d\n',i,nFiles);
    fprintf('%s\n',ftirFiles(i).name);
    fprintf('========================================\n');


    %% -----------------------------------------------------
    % Filename
    %% -----------------------------------------------------

    [~,baseName,~] = fileparts(ftirFiles(i).name);

    FileName(i) = string(baseName);

    baseLower = lower(string(baseName));


    %% -----------------------------------------------------
    % Cancer/non-cancer label
    %% -----------------------------------------------------

    if ismember(baseLower,nonCancerNames)

        Group(i) = "Non-cancer";

    elseif ismember(baseLower,cancerNames)

        Group(i) = "Cancer";

    else

        Group(i) = "Unknown";

    end


    %% -----------------------------------------------------
    % Locate H&E image
    %% -----------------------------------------------------

    heFile = fullfile(heFolder,[baseName '.png']);

    if ~exist(heFile,'file')

        warning('H&E image not found.');
        continue

    end


    %% =====================================================
    % STEP 1
    % READ FTIR
    %% =====================================================

    ftirFile = fullfile( ...
        ftirFolder, ...
        ftirFiles(i).name);

    cube = single(h5read(ftirFile,datasetName));


    if ~isequal(size(cube),[256 256 965])

        warning('Unexpected FTIR size.');
        continue

    end


    %% =====================================================
    % STEP 2
    % CREATE TISSUE MASK
    %% =====================================================

    Xall = reshape(cube,[],965);

    validPixels = ...
        all(isfinite(Xall),2) & ...
        any(abs(Xall) > 1e-12,2);

    tissueMask = reshape(validPixels,256,256);


    %% =====================================================
    % STEP 3
    % PCA64 SMOOTHING
    %% =====================================================

    X = double(Xall(validPixels,:));

    [coeff,score,~,~,~,mu] = ...
        pca(X,'NumComponents',FTIRsmoothPC);


    Xsmooth = ...
        score(:,1:FTIRsmoothPC) * ...
        coeff(:,1:FTIRsmoothPC)' + mu;


    %% =====================================================
    % STEP 4
    % FTIR SPECTRAL PREPROCESSING
    %% =====================================================

    Xprocessed = preprocessFTIR(Xsmooth);


    %% =====================================================
    % RESTORE FTIR CUBE
    %% =====================================================

    processedAll = ...
        zeros(256*256,965,'single');

    processedAll(validPixels,:) = ...
        single(Xprocessed);

    processedCube = ...
        reshape(processedAll,256,256,965);


    %% =====================================================
    % STEP 5
    % READ H&E
    %% =====================================================

    HE = imread(heFile);

    HE = im2single(HE);


    %% =====================================================
    % H&E 1024 x 1024 x 3
    %       ->
    %     256 x 256 x 3
    %% =====================================================

    HE256 = imresize(HE,[256 256],'bicubic');


    %% =====================================================
    % STEP 6
    % EXTRACT MATCHED 16x16 PATCHES
    %% =====================================================

    [Fpatch,Hpatch] = ...
        extractPatchFeatures( ...
        processedCube, ...
        HE256, ...
        tissueMask, ...
        patchSize, ...
        minimumTissueFraction);


    NumberPatches(i) = size(Fpatch,1);


    fprintf('Valid patches = %d\n', ...
        NumberPatches(i));


    %% -----------------------------------------------------
    % Need enough patches
    %% -----------------------------------------------------

    if size(Fpatch,1) < 30

        warning('Too few patches.');
        continue

    end


    %% =====================================================
    % STEP 7
    % PATCH FEATURE PCA
    %% =====================================================

    %
    % FTIR:
    % N x 965
    %       ->
    % N x 5
    %

    [~,Fscore] = ...
        pca(double(Fpatch), ...
        'NumComponents',patchFeaturePC);


    %
    % H&E:
    % N x 768
    %       ->
    % N x 5
    %

    [~,Hscore] = ...
        pca(double(Hpatch), ...
        'NumComponents',patchFeaturePC);


    F = Fscore(:,1:patchFeaturePC);

    H = Hscore(:,1:patchFeaturePC);


    %% =====================================================
    % STEP 8
    % MULTIVARIATE INFORMATION
    %% =====================================================

    metrics = ...
        multivariateInformation(F,H);


    MI(i) = metrics.MI;

    EntropyFTIR(i) = metrics.HF;

    EntropyHE(i) = metrics.HH;

    NormalizedMI(i) = metrics.NMI;

    SharedFTIR(i) = metrics.SharedFTIR;

    SharedHE(i) = metrics.SharedHE;


    fprintf('MI          = %.5f bits\n',MI(i));
    fprintf('NMI         = %.5f\n',NormalizedMI(i));
    fprintf('Shared FTIR = %.5f\n',SharedFTIR(i));
    fprintf('Shared H&E  = %.5f\n',SharedHE(i));

end


%% =========================================================
% RESULTS TABLE
%% =========================================================

Results = table( ...
    FileName, ...
    Group, ...
    NumberPatches, ...
    MI, ...
    EntropyFTIR, ...
    EntropyHE, ...
    NormalizedMI, ...
    SharedFTIR, ...
    SharedHE);


%% =========================================================
% SAVE
%% =========================================================




save( ...
    fullfile(outputFolder, ...
    'Patch_Multivariate_Information.mat'), ...
    'Results');




disp(Results);