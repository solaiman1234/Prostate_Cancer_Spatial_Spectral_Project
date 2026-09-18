clc
clear
close all

%% ============================================================
% INPUT FOLDERS
%% ============================================================

ftirFolder = ...
'D:\Prostate Cancer Assignment\spectral_annotation_965';

heFolder = ...
'D:\Prostate Cancer Assignment\aligned_he_annotation_only';


%% ============================================================
% OUTPUT FOLDER
%% ============================================================

outputFolder = ...
'D:\Prostate Cancer Assignment\MI_ConditionalEntropy_Results';

if ~exist(outputFolder,'dir')
    mkdir(outputFolder);
end


%% ============================================================
% Folder for FTIR energy maps
%% ============================================================

energyFolder = fullfile(outputFolder,'FTIR_Energy_Maps');

if ~exist(energyFolder,'dir')
    mkdir(energyFolder);
end


%% ============================================================
% Optional folder for resized H&E images
%% ============================================================

he256Folder = fullfile(outputFolder,'HE_256');

if ~exist(he256Folder,'dir')
    mkdir(he256Folder);
end


%% ============================================================
% OPTIONAL:
% Saving all reconstructed 256x256x965 cubes will require
% substantial disk space.
%
% false = reconstruct in memory but only save energy map
% true  = also save PCA64 reconstructed cube
%% ============================================================

saveReconstructedCube = false;

if saveReconstructedCube

    reconstructionFolder = ...
        fullfile(outputFolder,'PCA64_Reconstructed');

    if ~exist(reconstructionFolder,'dir')
        mkdir(reconstructionFolder);
    end

end


%% ============================================================
% PARAMETERS
%% ============================================================

numPC = 64;

numBins = 64;

nRows  = 256;
nCols  = 256;
nBands = 965;


%% ============================================================
% FIND ALL FTIR FILES
%% ============================================================

files = dir(fullfile(ftirFolder,'*.h5'));

nFiles = numel(files);

fprintf('Total FTIR files found = %d\n',nFiles);


%% ============================================================
% PREALLOCATE RESULTS
%% ============================================================

FileName = strings(nFiles,1);

ValidPixels = zeros(nFiles,1);

PCA64Variance = nan(nFiles,1);

Entropy_FTIR = nan(nFiles,1);

Entropy_HE = nan(nFiles,1);

JointEntropy = nan(nFiles,1);

MutualInformation = nan(nFiles,1);

NormalizedMI = nan(nFiles,1);

H_FTIR_given_HE = nan(nFiles,1);

H_HE_given_FTIR = nan(nFiles,1);

Normalized_H_FTIR_given_HE = nan(nFiles,1);

Normalized_H_HE_given_FTIR = nan(nFiles,1);


%% ============================================================
% PROCESS ALL FILES
%% ============================================================

for i = 1:nFiles

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('Processing %d / %d\n',i,nFiles);
    fprintf('%s\n',files(i).name);
    fprintf('============================================\n');

    FileName(i) = string(files(i).name);


    %% ========================================================
    % GET BASENAME
    %
    % example:
    %
    % ABC123.h5
    % ->
    % ABC123.png
    %% ========================================================

    [~,baseName,~] = fileparts(files(i).name);

    heFile = fullfile(heFolder,[baseName '.png']);


    if ~exist(heFile,'file')

        warning('H&E image not found for: %s',baseName);

        continue

    end


    %% ========================================================
    % LOAD FTIR
    %% ========================================================

    ftirFile = fullfile(ftirFolder,files(i).name);

    FTIR = h5read(ftirFile,'/spectra');

    FTIR = single(FTIR);


    %% ========================================================
    % CHECK SIZE
    %% ========================================================

    if ~isequal(size(FTIR),[256 256 965])

        warning('Unexpected FTIR size for %s', ...
            files(i).name);

        disp(size(FTIR));

        continue

    end


    %% ========================================================
    % RESHAPE
    %
    % 256 x 256 x 965
    %
    % ->
    %
    % 65536 x 965
    %% ========================================================

    Xall = reshape(FTIR,[],nBands);


    %% ========================================================
    % IDENTIFY VALID FTIR TISSUE PIXELS
    %
    % Background assumed to be zero.
    %% ========================================================

    validMask = ...
        all(isfinite(Xall),2) & ...
        any(abs(Xall) > 0,2);


    X = Xall(validMask,:);


    ValidPixels(i) = size(X,1);


    fprintf('Valid tissue pixels = %d\n', ...
        ValidPixels(i));


    if size(X,1) <= numPC

        warning('Not enough valid pixels for PCA.');

        continue

    end


    %% ========================================================
    % PCA 64 SMOOTHING
    %% ========================================================

    fprintf('Running PCA64...\n');


    [coeff,score,latent,~,explained,mu] = ...
        pca(X, ...
            'NumComponents',numPC, ...
            'Algorithm','eig');


    PCA64Variance(i) = sum(explained);


    fprintf('Variance represented by PCA64 = %.3f %%\n', ...
        PCA64Variance(i));


    %% ========================================================
    % RECONSTRUCT BACK TO 965 CHANNELS
    %
    % N x 64
    %   x
    % 64 x 965
    %
    % ->
    %
    % N x 965
    %% ========================================================

    X_reconstructed = ...
        score(:,1:numPC) * ...
        coeff(:,1:numPC)' + mu;


    %% ========================================================
    % CREATE FTIR ENERGY MAP
    %
    % RMS spectral energy:
    %
    % sqrt(mean(spectrum.^2))
    %
    % This produces one scalar value for every spatial pixel.
    %% ========================================================

    energyTissue = sqrt( ...
        mean(X_reconstructed.^2,2));


    %% ========================================================
    % PUT TISSUE VALUES BACK INTO 256x256 IMAGE
    %% ========================================================

    energyVector = ...
        zeros(nRows*nCols,1,'single');


    energyVector(validMask) = ...
        single(energyTissue);


    energyMap = reshape( ...
        energyVector, ...
        nRows,nCols);


    %% ========================================================
    % SAVE ENERGY MAP
    %
    % Use same basename
    %% ========================================================

    save( ...
        fullfile(energyFolder,[baseName '.mat']), ...
        'energyMap');


    %% ========================================================
    % OPTIONAL:
    % Save PCA64 reconstructed hyperspectral cube
    %% ========================================================

    if saveReconstructedCube

        reconstructedAll = ...
            zeros(nRows*nCols,nBands,'single');

        reconstructedAll(validMask,:) = ...
            single(X_reconstructed);

        reconstructedCube = reshape( ...
            reconstructedAll, ...
            nRows,nCols,nBands);


        outputH5 = fullfile( ...
            reconstructionFolder, ...
            files(i).name);


        if exist(outputH5,'file')
            delete(outputH5);
        end


        h5create( ...
            outputH5, ...
            '/spectra', ...
            size(reconstructedCube), ...
            'Datatype','single', ...
            'ChunkSize',[64 64 32], ...
            'Deflate',4);


        h5write( ...
            outputH5, ...
            '/spectra', ...
            reconstructedCube);

    end


    %% ========================================================
    % LOAD PAIRED H&E IMAGE
    %% ========================================================

    HE = imread(heFile);


    %% ========================================================
    % HANDLE RGB OR GRAYSCALE
    %% ========================================================

    if size(HE,3) == 3

        HE = im2single(HE);

        HE256_RGB = imresize( ...
            HE,[256 256],'bicubic');

        HE256 = rgb2gray(HE256_RGB);

    else

        HE = im2single(HE);

        HE256 = imresize( ...
            HE,[256 256],'bicubic');

    end


    %% ========================================================
    % SAVE RESIZED H&E FOR VISUAL CHECK
    %% ========================================================

    imwrite( ...
        mat2gray(HE256), ...
        fullfile(he256Folder,[baseName '.png']));


    %% ========================================================
    % MAKE SPATIAL FTIR MASK
    %% ========================================================

    validMask2D = reshape( ...
        validMask,nRows,nCols);


    %% ========================================================
    % USE EXACTLY THE SAME PIXELS FROM BOTH MODALITIES
    %% ========================================================

    F = double(energyMap(validMask2D));

    H = double(HE256(validMask2D));


    %% ========================================================
    % REMOVE ANY INVALID VALUES
    %% ========================================================

    good = ...
        isfinite(F) & ...
        isfinite(H);


    F = F(good);

    H = H(good);


    %% ========================================================
    % ROBUST NORMALIZATION
    %
    % Use 1st and 99th percentiles to reduce outlier effects.
    %% ========================================================

    F = robustNormalize01(F);

    H = robustNormalize01(H);


    %% ========================================================
    % MUTUAL INFORMATION + CONDITIONAL ENTROPY
    %% ========================================================

    metrics = calculateInformationMetrics( ...
        F,H,numBins);


    %% ========================================================
    % SAVE METRICS
    %% ========================================================

    Entropy_FTIR(i) = metrics.HF;

    Entropy_HE(i) = metrics.HH;

    JointEntropy(i) = metrics.HFH;

    MutualInformation(i) = metrics.MI;

    NormalizedMI(i) = metrics.NMI;

    H_FTIR_given_HE(i) = ...
        metrics.H_F_given_H;

    H_HE_given_FTIR(i) = ...
        metrics.H_H_given_F;

    Normalized_H_FTIR_given_HE(i) = ...
        metrics.Normalized_H_F_given_H;

    Normalized_H_HE_given_FTIR(i) = ...
        metrics.Normalized_H_H_given_F;


    %% ========================================================
    % DISPLAY CURRENT RESULT
    %% ========================================================

    fprintf('\n');

    fprintf('MI = %.4f bits\n', ...
        MutualInformation(i));

    fprintf('NMI = %.4f\n', ...
        NormalizedMI(i));

    fprintf('H(FTIR | H&E) = %.4f bits\n', ...
        H_FTIR_given_HE(i));

    fprintf('H(H&E | FTIR) = %.4f bits\n', ...
        H_HE_given_FTIR(i));


    %% ========================================================
    % CLEAR LARGE VARIABLES BEFORE NEXT FILE
    %% ========================================================

    clear FTIR Xall X coeff score latent mu
    clear X_reconstructed

end


%% ============================================================
% CREATE RESULTS TABLE
%% ============================================================

Results = table( ...
    FileName, ...
    ValidPixels, ...
    PCA64Variance, ...
    Entropy_FTIR, ...
    Entropy_HE, ...
    JointEntropy, ...
    MutualInformation, ...
    NormalizedMI, ...
    H_FTIR_given_HE, ...
    H_HE_given_FTIR, ...
    Normalized_H_FTIR_given_HE, ...
    Normalized_H_HE_given_FTIR);







%% ============================================================
% SAVE MATLAB DATA
%% ============================================================

save( ...
    fullfile(outputFolder, ...
    'Multimodal_MI_ConditionalEntropy.mat'), ...
    'Results');


fprintf('\n');
fprintf('============================================\n');
fprintf('ALL FILES COMPLETED\n');
fprintf('============================================\n');

fprintf('Results saved to:\n%s\n',csvFile);