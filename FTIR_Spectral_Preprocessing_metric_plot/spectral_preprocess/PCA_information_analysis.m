clc
clear
close all


%% =========================================================
% Input folder
%% =========================================================

folderPath = ...
'D:\Prostate Cancer Assignment\spectral_annotation_965';


%% =========================================================
% Output folder
%% =========================================================

outputFolder = ...
'D:\Prostate Cancer Assignment\PCA_information_analysis';

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end


%% =========================================================
% Find all H5 files
%% =========================================================

files = dir(fullfile(folderPath, '*.h5'));

nFiles = numel(files);

fprintf('Total number of files = %d\n', nFiles);


%% =========================================================
% PCA levels
%
% First entry means ALL 965 channels
%% =========================================================

K_values = [965 256 64 16 4 1];

K_names = ["All", "256", "64", "16", "4", "1"];

nLevels = length(K_values);


%% =========================================================
% Savitzky-Golay parameters
%% =========================================================

polyOrder = 5;
frameLength = 17;

halfWindow = (frameLength - 1) / 2;


%% =========================================================
% Data dimensions
%% =========================================================

nRows = 256;
nCols = 256;
nBands = 965;

validChannels = ...
    (halfWindow + 1):(nBands - halfWindow);


%% =========================================================
% Preallocate results
%% =========================================================

fileNames = strings(nFiles,1);

validPixelCount = zeros(nFiles,1);


% PCA variance retained
varianceRetained = nan(nFiles,nLevels);


% Reconstruction metrics
NRMSE = nan(nFiles,nLevels);

meanSAM = nan(nFiles,nLevels);

medianSAM = nan(nFiles,nLevels);


% Final processed spectrum comparison
derivativeCorrelation = nan(nFiles,nLevels);


% Spatial preservation
energySSIM = nan(nFiles,nLevels);

edgeDice = nan(nFiles,nLevels);


%% =========================================================
% Mean processed spectrum from every image
%
% Dimensions:
%
% 227 × 6 × 965
%
% This will be very useful later when you separate
% cancer and normal cases.
%% =========================================================

meanProcessedSpectrum = ...
    nan(nFiles,nLevels,nBands,'single');


%% =========================================================
% MAIN LOOP
%% =========================================================

for i = 1:nFiles

    fprintf('\n');
    fprintf('=========================================\n');
    fprintf('Processing %d of %d\n', i, nFiles);
    fprintf('%s\n', files(i).name);
    fprintf('=========================================\n');


    fileNames(i) = string(files(i).name);


    %% -----------------------------------------------------
    % Load H5
    %% -----------------------------------------------------

    fullFileName = ...
        fullfile(folderPath, files(i).name);

    cube = h5read(fullFileName, '/spectra');

    cube = single(cube);


    %% -----------------------------------------------------
    % Check dimension
    %% -----------------------------------------------------

    if ~isequal(size(cube), [256 256 965])

        warning('Unexpected size in %s', files(i).name);

        continue;

    end


    %% -----------------------------------------------------
    % Reshape:
    %
    % 256 × 256 × 965
    %
    % becomes
    %
    % 65536 × 965
    %% -----------------------------------------------------

    Xall = reshape(cube, [], nBands);


    %% =====================================================
    % Remove invalid / background pixels
    %
    % Assumes annotation background is zero.
    %% =====================================================

    finitePixels = all(isfinite(Xall), 2);

    nonZeroPixels = ...
        sum(abs(Xall), 2) > 0;

    validMask = ...
        finitePixels & nonZeroPixels;


    X = Xall(validMask, :);


    validPixelCount(i) = size(X,1);

    fprintf('Valid pixels = %d\n', size(X,1));


    if size(X,1) <= 256

        warning('Too few valid pixels in %s', ...
                files(i).name);

        continue;

    end


    %% =====================================================
    % ALL COMPONENT REFERENCE
    %
    % No PCA
    %
    % Raw spectrum
    %     ↓
    % Shift
    %     ↓
    % Normalize
    %     ↓
    % SG second derivative
    %% =====================================================

    fprintf('Processing ALL components...\n');

    D_reference = ...
        preprocess_FTIR(X, ...
                        polyOrder, ...
                        frameLength);


    %% -----------------------------------------------------
    % Save mean processed spectrum
    %% -----------------------------------------------------

    meanProcessedSpectrum(i,1,:) = ...
        mean(D_reference, 1, 'omitnan');


    %% -----------------------------------------------------
    % ALL = perfect reference
    %% -----------------------------------------------------

    varianceRetained(i,1) = 100;

    NRMSE(i,1) = 0;

    meanSAM(i,1) = 0;

    medianSAM(i,1) = 0;

    derivativeCorrelation(i,1) = 1;

    energySSIM(i,1) = 1;

    edgeDice(i,1) = 1;


    %% =====================================================
    % Create reference spectral-energy image
    %% =====================================================

    energyReference = sqrt( ...
        sum(D_reference(:,validChannels).^2, 2));


    energyVector = ...
        zeros(nRows*nCols, 1, 'single');


    energyVector(validMask) = ...
        energyReference;


    energyMapReference = ...
        reshape(energyVector, nRows, nCols);


    %% =====================================================
    % Tissue bounding box
    %% =====================================================

    tissueMask2D = ...
        reshape(validMask, nRows, nCols);


    [r,c] = find(tissueMask2D);


    if isempty(r)

        warning('No tissue found');

        continue;

    end


    row1 = min(r);
    row2 = max(r);

    col1 = min(c);
    col2 = max(c);


    referenceCrop = ...
        energyMapReference(row1:row2, ...
                           col1:col2);


    %% =====================================================
    % PCA
    %
    % ONLY CALCULATE ONCE.
    %
    % Largest reduced representation = 256 PCs
    %% =====================================================

    fprintf('Calculating PCA...\n');


    [coeff, score, latent, ~, ~, mu] = ...
        pca(X, ...
            'NumComponents',256, ...
            'Algorithm','eig');


    %% -----------------------------------------------------
    % Total variance in original data
    %% -----------------------------------------------------

    totalVariance = ...
        sum(var(double(X), 0, 1));


    %% =====================================================
    % LOOP THROUGH
    %
    % 256
    % 64
    % 16
    % 4
    % 1
    %% =====================================================

    for level = 2:nLevels

        K = K_values(level);


        fprintf('   Reconstructing PCA %d...\n', K);


        %% =============================================
        % PCA reconstruction
        %% =============================================

        X_reconstructed = ...
            score(:,1:K) * ...
            coeff(:,1:K)' + mu;


        %% =============================================
        % Variance retained
        %% =============================================

        varianceRetained(i,level) = ...
            100 * sum(latent(1:K)) / ...
            totalVariance;


        %% =============================================
        % NRMSE
        %
        % ORIGINAL vs PCA reconstruction
        %% =============================================

        errorData = ...
            X - X_reconstructed;


        RMSE_value = ...
            sqrt(mean(errorData(:).^2));


        signalRange = ...
            max(X(:)) - min(X(:));


        NRMSE(i,level) = ...
            RMSE_value / ...
            (signalRange + eps);


        %% =============================================
        % Spectral Angle Mapper
        %% =============================================

        dotValue = ...
            sum(X .* X_reconstructed, 2);


        originalNorm = ...
            sqrt(sum(X.^2,2));


        reconstructedNorm = ...
            sqrt(sum(X_reconstructed.^2,2));


        cosTheta = ...
            dotValue ./ ...
            (originalNorm .* ...
             reconstructedNorm + eps);


        cosTheta = ...
            max(min(cosTheta,1),-1);


        SAM = acosd(cosTheta);


        meanSAM(i,level) = ...
            mean(SAM,'omitnan');


        medianSAM(i,level) = ...
            median(SAM,'omitnan');


        %% =============================================
        % Apply preprocessing
        %
        % PCA reconstruction
        %       ↓
        % Shift
        %       ↓
        % Normalize
        %       ↓
        % SG 2nd derivative
        %% =============================================

        D_PCA = ...
            preprocess_FTIR( ...
                X_reconstructed, ...
                polyOrder, ...
                frameLength);


        %% =============================================
        % Save mean processed spectrum
        %% =============================================

        meanProcessedSpectrum(i,level,:) = ...
            mean(D_PCA,1,'omitnan');


        %% =============================================
        % Correlation with ALL-component
        % processed spectrum
        %% =============================================

        R_pixel = ...
            rowCorrelation( ...
                D_reference(:,validChannels), ...
                D_PCA(:,validChannels));


        derivativeCorrelation(i,level) = ...
            mean(R_pixel,'omitnan');


        %% =============================================
        % Spectral energy map
        %% =============================================

        energyPCA = sqrt( ...
            sum(D_PCA(:,validChannels).^2,2));


        energyVectorPCA = ...
            zeros(nRows*nCols,1,'single');


        energyVectorPCA(validMask) = ...
            energyPCA;


        energyMapPCA = ...
            reshape(energyVectorPCA, ...
                    nRows,nCols);


        %% =============================================
        % Crop to tissue region
        %% =============================================

        PCAcrop = ...
            energyMapPCA(row1:row2, ...
                         col1:col2);


        %% =============================================
        % SSIM
        %% =============================================

        dynamicRange = ...
            max(referenceCrop(:)) - ...
            min(referenceCrop(:));


        if dynamicRange > 0

            energySSIM(i,level) = ...
                ssim( ...
                    PCAcrop, ...
                    referenceCrop, ...
                    'DynamicRange', ...
                    double(dynamicRange));

        end


        %% =============================================
        % Edge Dice
        %% =============================================

        referenceEdge = ...
            edge(mat2gray(referenceCrop), ...
                 'Canny');


        PCAedge = ...
            edge(mat2gray(PCAcrop), ...
                 'Canny');


        denominator = ...
            sum(referenceEdge(:)) + ...
            sum(PCAedge(:));


        if denominator > 0

            edgeDice(i,level) = ...
                2 * sum( ...
                    referenceEdge(:) & ...
                    PCAedge(:)) ...
                / denominator;

        end


        %% =============================================
        % Clear large arrays
        %% =============================================

        clear X_reconstructed
        clear D_PCA
        clear energyPCA
        clear energyVectorPCA
        clear energyMapPCA

    end


    %% =====================================================
    % Clear image before next file
    %% =====================================================

    clear cube
    clear Xall
    clear X
    clear coeff
    clear score
    clear latent
    clear mu
    clear D_reference

end


%% =========================================================
% SAVE EVERYTHING
%% =========================================================

save(fullfile(outputFolder, ...
    'PCA_information_results.mat'), ...
    'fileNames', ...
    'K_values', ...
    'K_names', ...
    'validPixelCount', ...
    'varianceRetained', ...
    'NRMSE', ...
    'meanSAM', ...
    'medianSAM', ...
    'derivativeCorrelation', ...
    'energySSIM', ...
    'edgeDice', ...
    'meanProcessedSpectrum', ...
    'validChannels', ...
    '-v7.3');


fprintf('\n');
fprintf('=========================================\n');
fprintf('ALL 227 FILES FINISHED\n');
fprintf('=========================================\n');