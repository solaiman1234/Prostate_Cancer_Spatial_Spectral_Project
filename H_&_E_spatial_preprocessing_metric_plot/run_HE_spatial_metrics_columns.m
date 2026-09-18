function Results = run_HE_spatial_metrics_columns()

clc;
close all;

%% =========================================================
% H&E SPATIAL INFORMATION ANALYSIS
%
% Rows    = individual H&E images
% Columns = [1024 512 256 128 64]
%
% The lower-resolution images are assumed to have already
% been upsampled back to 1024x1024 for comparison.
%
% Background is excluded using the corresponding 256x256
% annotation mask, resized to 1024x1024.
%% =========================================================


%% =========================================================
% 1. FOLDERS
%% =========================================================

rootFolder = ...
    'D:\Prostate Cancer Assignment\HE_spatial_resolution';

referenceFolder = ...
    fullfile(rootFolder,'Resolution_1024');

comparisonRoot = ...
    fullfile(rootFolder,'Upsampled_to_1024');

maskFolder = ...
    'D:\Prostate Cancer Assignment\mask_annotation_only';

resultFolder = ...
    fullfile(rootFolder,'Metric_Results_Tissue_Only');


if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end


%% =========================================================
% 2. RESOLUTION ORDER
%
% IMPORTANT:
% Every metric matrix follows this column order.
%% =========================================================

resolutionList = [1024 512 256 128 64];

nResolution = length(resolutionList);


%% =========================================================
% 3. FIND REFERENCE H&E IMAGES
%% =========================================================

filesPNG  = dir(fullfile(referenceFolder,'*.png'));


files = [ ...
    filesPNG;
    ];


nFiles = length(files);

fprintf('Total H&E images found = %d\n',nFiles);


if nFiles == 0
    error('No H&E images were found in the reference folder.');
end


%% =========================================================
% 4. FILE NAMES
%% =========================================================

FileName = cell(nFiles,1);

MaskFileName = cell(nFiles,1);


%% =========================================================
% 5. PREALLOCATE METRICS
%
% N images x 5 spatial resolutions
%% =========================================================

SSIM_value = nan(nFiles,nResolution);

PSNR_value = nan(nFiles,nResolution);

RMSE_value = nan(nFiles,nResolution);

NRMSE_value = nan(nFiles,nResolution);

MAE_value = nan(nFiles,nResolution);

ImageCorrelation = nan(nFiles,nResolution);

EdgeDice = nan(nFiles,nResolution);


%% Fine morphology

GradientEnergy = nan(nFiles,nResolution);

GradientRetention = nan(nFiles,nResolution);

LaplacianVariance = nan(nFiles,nResolution);

LaplacianRetention = nan(nFiles,nResolution);


%% Complexity

EntropyValue = nan(nFiles,nResolution);

EntropyRetention = nan(nFiles,nResolution);


%% Spatial-frequency information

HFPower = nan(nFiles,nResolution);

HFRetention = nan(nFiles,nResolution);


%% Texture

GLCMContrast = nan(nFiles,nResolution);

GLCMCorrelation = nan(nFiles,nResolution);

GLCMEnergy = nan(nFiles,nResolution);

GLCMHomogeneity = nan(nFiles,nResolution);


%% Tissue information

TissuePixelCount = nan(nFiles,1);

TissueFraction = nan(nFiles,1);


%% =========================================================
% 6. PROCESS EVERY IMAGE
%% =========================================================

for i = 1:nFiles

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('Processing image %d / %d\n',i,nFiles);
    fprintf('%s\n',files(i).name);
    fprintf('============================================\n');


    FileName{i} = files(i).name;


    %% =====================================================
    % LOAD 1024 REFERENCE IMAGE
    %% =====================================================

    referenceFile = ...
        fullfile(referenceFolder,files(i).name);

    Iref = imread(referenceFile);

    Iref = makeRGB(Iref);

    Iref = im2double(Iref);


    if size(Iref,1) ~= 1024 || size(Iref,2) ~= 1024

        warning('Reference image is not 1024x1024: %s', ...
            files(i).name);

        continue;

    end


    %% =====================================================
    % FIND CORRESPONDING MASK
    %% =====================================================

    maskFile = ...
        findMaskFile(maskFolder,files(i).name);


    if isempty(maskFile)

        warning('Mask not found for %s',files(i).name);

        continue;

    end


    [~,maskName,maskExt] = fileparts(maskFile);

    MaskFileName{i} = [maskName maskExt];


    %% =====================================================
    % LOAD MASK
    %% =====================================================

    M = imread(maskFile);


    if ndims(M) == 3
        M = M(:,:,1);
    end


    %% -----------------------------------------------------
    % Assumption:
    %
    % background = 0
    % tissue     = non-zero
    %% -----------------------------------------------------

    M = M > 0;


    %% =====================================================
    % RESIZE MASK:
    %
    % 256x256 -> 1024x1024
    %
    % nearest neighbour MUST be used for categorical masks.
    %% =====================================================

    tissueMask = ...
        imresize(M,[1024 1024],'nearest');

    tissueMask = logical(tissueMask);


    %% =====================================================
    % ERODE SLIGHTLY
    %
    % Removes the artificial tissue/background boundary
    % from edge/gradient calculations.
    %% =====================================================

    analysisMask = ...
        imerode(tissueMask,ones(5,5));


    %% Safety check

    if nnz(analysisMask) < 0.5*nnz(tissueMask)

        analysisMask = tissueMask;

    end


    TissuePixelCount(i) = ...
        nnz(analysisMask);


    TissueFraction(i) = ...
        nnz(analysisMask) / numel(analysisMask);


    if TissuePixelCount(i) < 100

        warning('Too few valid tissue pixels: %s', ...
            files(i).name);

        continue;

    end


    %% =====================================================
    % REFERENCE GRAYSCALE
    %% =====================================================

    grayRef = rgb2gray(Iref);


    %% =====================================================
    % 7. CALCULATE REFERENCE FEATURES ONCE
    %% =====================================================


    %% -----------------------------------------------------
    % REFERENCE GRADIENT ENERGY
    %% -----------------------------------------------------

    [GxRef,GyRef] = ...
        imgradientxy(grayRef,'sobel');


    gradientMapRef = ...
        GxRef.^2 + GyRef.^2;


    referenceGradientPixels = ...
        gradientMapRef(analysisMask);


    gradientRef = ...
        mean(referenceGradientPixels);


    %% -----------------------------------------------------
    % REFERENCE LAPLACIAN VARIANCE
    %% -----------------------------------------------------

    lapFilter = ...
        fspecial('laplacian',0.2);


    lapRef = ...
        imfilter(grayRef,lapFilter,'replicate');


    lapPixelsRef = ...
        lapRef(analysisMask);


    laplacianRef = ...
        var(lapPixelsRef);


    %% -----------------------------------------------------
    % REFERENCE ENTROPY
    %% -----------------------------------------------------

    entropyRef = ...
        calculateMaskedEntropy( ...
        grayRef,analysisMask);


    %% -----------------------------------------------------
    % REFERENCE HIGH-FREQUENCY POWER
    %% -----------------------------------------------------

    HFref = ...
        calculateMaskedHFPower( ...
        grayRef,analysisMask);


    %% -----------------------------------------------------
    % REFERENCE GLCM
    %% -----------------------------------------------------

    glcmRef = ...
        calculateMaskedGLCM( ...
        grayRef,analysisMask,16);


    %% -----------------------------------------------------
    % REFERENCE EDGE MAP
    %% -----------------------------------------------------

    [edgeRef,cannyThreshold] = ...
        edge(grayRef,'Canny');


    edgeRef = ...
        edgeRef & analysisMask;


    %% =====================================================
    % 8. PROCESS EACH SPATIAL RESOLUTION
    %% =====================================================

    for r = 1:nResolution

        currentResolution = ...
            resolutionList(r);


        fprintf('   Resolution %d\n',currentResolution);


        %% =================================================
        % LOAD COMPARISON IMAGE
        %% =================================================

        if currentResolution == 1024

            Icomp = Iref;

        else

            comparisonFolder = ...
                fullfile( ...
                comparisonRoot, ...
                ['From_' num2str(currentResolution)]);


            comparisonFile = ...
                fullfile( ...
                comparisonFolder, ...
                files(i).name);


            if ~exist(comparisonFile,'file')

                warning('Comparison image not found: %s', ...
                    comparisonFile);

                continue;

            end


            Icomp = imread(comparisonFile);

            Icomp = makeRGB(Icomp);

            Icomp = im2double(Icomp);

        end


        %% -------------------------------------------------
        % Verify comparison image size
        %% -------------------------------------------------

        if size(Icomp,1) ~= 1024 || ...
           size(Icomp,2) ~= 1024

            warning('Comparison image is not 1024x1024.');

            continue;

        end


        %% =================================================
        % GRAYSCALE
        %% =================================================

        grayComp = ...
            rgb2gray(Icomp);


        %% =================================================
        % TISSUE-ONLY PIXELS
        %% =================================================

        refPixels = ...
            grayRef(analysisMask);


        compPixels = ...
            grayComp(analysisMask);


        differencePixels = ...
            refPixels - compPixels;


        %% =================================================
        % RMSE
        %% =================================================

        mseValue = ...
            mean(differencePixels.^2);


        RMSE_value(i,r) = ...
            sqrt(mseValue);


        %% =================================================
        % NRMSE
        %% =================================================

        referenceRange = ...
            max(refPixels) - min(refPixels);


        if referenceRange > 0

            NRMSE_value(i,r) = ...
                RMSE_value(i,r) / referenceRange;

        else

            NRMSE_value(i,r) = NaN;

        end


        %% =================================================
        % MAE
        %% =================================================

        MAE_value(i,r) = ...
            mean(abs(differencePixels));


        %% =================================================
        % PSNR
        %% =================================================

        if mseValue == 0

            PSNR_value(i,r) = Inf;

        else

            PSNR_value(i,r) = ...
                10*log10(1/mseValue);

        end


        %% =================================================
        % IMAGE CORRELATION
        %% =================================================

        if std(refPixels) > 0 && ...
           std(compPixels) > 0

            C = ...
                corrcoef(refPixels,compPixels);


            ImageCorrelation(i,r) = ...
                C(1,2);

        else

            ImageCorrelation(i,r) = NaN;

        end


        %% =================================================
        % MASKED SSIM
        %% =================================================

        SSIM_value(i,r) = ...
            calculateMaskedSSIM( ...
            grayRef, ...
            grayComp, ...
            analysisMask);


        %% =================================================
        % EDGE DICE
        %% =================================================

        edgeComp = ...
            edge( ...
            grayComp, ...
            'Canny', ...
            cannyThreshold);


        edgeComp = ...
            edgeComp & analysisMask;


        edgeDenominator = ...
            nnz(edgeRef) + nnz(edgeComp);


        if edgeDenominator > 0

            EdgeDice(i,r) = ...
                2*nnz(edgeRef & edgeComp) / ...
                edgeDenominator;

        else

            EdgeDice(i,r) = NaN;

        end


        %% =================================================
        % GRADIENT ENERGY
        %% =================================================

        [GxComp,GyComp] = ...
            imgradientxy(grayComp,'sobel');


        gradientMapComp = ...
            GxComp.^2 + GyComp.^2;


        comparisonGradientPixels = ...
            gradientMapComp(analysisMask);


        gradientComp = ...
            mean(comparisonGradientPixels);


        GradientEnergy(i,r) = ...
            gradientComp;


        if gradientRef > 0

            GradientRetention(i,r) = ...
                gradientComp / gradientRef;

        else

            GradientRetention(i,r) = NaN;

        end


        %% =================================================
        % LAPLACIAN VARIANCE
        %% =================================================

        lapComp = ...
            imfilter( ...
            grayComp, ...
            lapFilter, ...
            'replicate');


        lapPixelsComp = ...
            lapComp(analysisMask);


        laplacianComp = ...
            var(lapPixelsComp);


        LaplacianVariance(i,r) = ...
            laplacianComp;


        if laplacianRef > 0

            LaplacianRetention(i,r) = ...
                laplacianComp / laplacianRef;

        else

            LaplacianRetention(i,r) = NaN;

        end


        %% =================================================
        % ENTROPY
        %% =================================================

        entropyComp = ...
            calculateMaskedEntropy( ...
            grayComp,analysisMask);


        EntropyValue(i,r) = ...
            entropyComp;


        if entropyRef > 0

            EntropyRetention(i,r) = ...
                entropyComp / entropyRef;

        else

            EntropyRetention(i,r) = NaN;

        end


        %% =================================================
        % HIGH-FREQUENCY FOURIER POWER
        %% =================================================

        HFcomp = ...
            calculateMaskedHFPower( ...
            grayComp,analysisMask);


        HFPower(i,r) = ...
            HFcomp;


        if HFref > 0

            HFRetention(i,r) = ...
                HFcomp / HFref;

        else

            HFRetention(i,r) = NaN;

        end


        %% =================================================
        % GLCM
        %% =================================================

        glcmComp = ...
            calculateMaskedGLCM( ...
            grayComp, ...
            analysisMask, ...
            16);


        GLCMContrast(i,r) = ...
            glcmComp.Contrast;


        GLCMCorrelation(i,r) = ...
            glcmComp.Correlation;


        GLCMEnergy(i,r) = ...
            glcmComp.Energy;


        GLCMHomogeneity(i,r) = ...
            glcmComp.Homogeneity;

    end

end


%% =========================================================
% 9. COLUMN LABELS
%% =========================================================

ColumnNames = { ...
    'FileName', ...
    'Resolution_1024', ...
    'Resolution_512', ...
    'Resolution_256', ...
    'Resolution_128', ...
    'Resolution_64'};


%% =========================================================
% 10. CREATE N x 6 CELL ARRAYS
%
% Column 1 = filename
% Column 2 = 1024
% Column 3 = 512
% Column 4 = 256
% Column 5 = 128
% Column 6 = 64
%% =========================================================

SSIM_data = ...
    [FileName num2cell(SSIM_value)];


PSNR_data = ...
    [FileName num2cell(PSNR_value)];


RMSE_data = ...
    [FileName num2cell(RMSE_value)];


NRMSE_data = ...
    [FileName num2cell(NRMSE_value)];


MAE_data = ...
    [FileName num2cell(MAE_value)];


ImageCorrelation_data = ...
    [FileName num2cell(ImageCorrelation)];


EdgeDice_data = ...
    [FileName num2cell(EdgeDice)];


GradientEnergy_data = ...
    [FileName num2cell(GradientEnergy)];


GradientRetention_data = ...
    [FileName num2cell(GradientRetention)];


LaplacianVariance_data = ...
    [FileName num2cell(LaplacianVariance)];


LaplacianRetention_data = ...
    [FileName num2cell(LaplacianRetention)];


Entropy_data = ...
    [FileName num2cell(EntropyValue)];


EntropyRetention_data = ...
    [FileName num2cell(EntropyRetention)];


HFPower_data = ...
    [FileName num2cell(HFPower)];


HFRetention_data = ...
    [FileName num2cell(HFRetention)];


GLCMContrast_data = ...
    [FileName num2cell(GLCMContrast)];


GLCMCorrelation_data = ...
    [FileName num2cell(GLCMCorrelation)];


GLCMEnergy_data = ...
    [FileName num2cell(GLCMEnergy)];


GLCMHomogeneity_data = ...
    [FileName num2cell(GLCMHomogeneity)];


%% =========================================================
% 11. CREATE RESULTS STRUCTURE
%% =========================================================

Results = struct;


Results.FileName = ...
    FileName;


Results.MaskFileName = ...
    MaskFileName;


Results.Resolution = ...
    resolutionList;


Results.ColumnNames = ...
    ColumnNames;


%% Conventional metrics

Results.SSIM = ...
    SSIM_value;

Results.PSNR = ...
    PSNR_value;

Results.RMSE = ...
    RMSE_value;

Results.NRMSE = ...
    NRMSE_value;

Results.MAE = ...
    MAE_value;

Results.ImageCorrelation = ...
    ImageCorrelation;


%% Morphology

Results.EdgeDice = ...
    EdgeDice;

Results.GradientEnergy = ...
    GradientEnergy;

Results.GradientRetention = ...
    GradientRetention;

Results.LaplacianVariance = ...
    LaplacianVariance;

Results.LaplacianRetention = ...
    LaplacianRetention;


%% Complexity

Results.Entropy = ...
    EntropyValue;

Results.EntropyRetention = ...
    EntropyRetention;


%% Frequency

Results.HFPower = ...
    HFPower;

Results.HFRetention = ...
    HFRetention;


%% Texture

Results.GLCMContrast = ...
    GLCMContrast;

Results.GLCMCorrelation = ...
    GLCMCorrelation;

Results.GLCMEnergy = ...
    GLCMEnergy;

Results.GLCMHomogeneity = ...
    GLCMHomogeneity;


%% Tissue

Results.TissuePixelCount = ...
    TissuePixelCount;

Results.TissueFraction = ...
    TissueFraction;


%% Combined filename + resolution values

Results.SSIM_data = ...
    SSIM_data;

Results.PSNR_data = ...
    PSNR_data;

Results.RMSE_data = ...
    RMSE_data;

Results.NRMSE_data = ...
    NRMSE_data;

Results.MAE_data = ...
    MAE_data;

Results.ImageCorrelation_data = ...
    ImageCorrelation_data;

Results.EdgeDice_data = ...
    EdgeDice_data;

Results.GradientEnergy_data = ...
    GradientEnergy_data;

Results.GradientRetention_data = ...
    GradientRetention_data;

Results.LaplacianVariance_data = ...
    LaplacianVariance_data;

Results.LaplacianRetention_data = ...
    LaplacianRetention_data;

Results.Entropy_data = ...
    Entropy_data;

Results.EntropyRetention_data = ...
    EntropyRetention_data;

Results.HFPower_data = ...
    HFPower_data;

Results.HFRetention_data = ...
    HFRetention_data;

Results.GLCMContrast_data = ...
    GLCMContrast_data;

Results.GLCMCorrelation_data = ...
    GLCMCorrelation_data;

Results.GLCMEnergy_data = ...
    GLCMEnergy_data;

Results.GLCMHomogeneity_data = ...
    GLCMHomogeneity_data;


%% =========================================================
% 12. SAVE MAT FILE
%% =========================================================

outputFile = fullfile( ...
    resultFolder, ...
    'HE_spatial_metrics_by_resolution.mat');


save( ...
    outputFile, ...
    'Results', ...
    'FileName', ...
    'MaskFileName', ...
    'resolutionList', ...
    'ColumnNames', ...
    'SSIM_value', ...
    'PSNR_value', ...
    'RMSE_value', ...
    'NRMSE_value', ...
    'MAE_value', ...
    'ImageCorrelation', ...
    'EdgeDice', ...
    'GradientEnergy', ...
    'GradientRetention', ...
    'LaplacianVariance', ...
    'LaplacianRetention', ...
    'EntropyValue', ...
    'EntropyRetention', ...
    'HFPower', ...
    'HFRetention', ...
    'GLCMContrast', ...
    'GLCMCorrelation', ...
    'GLCMEnergy', ...
    'GLCMHomogeneity', ...
    'TissuePixelCount', ...
    'TissueFraction', ...
    'SSIM_data', ...
    'PSNR_data', ...
    'RMSE_data', ...
    'NRMSE_data', ...
    'MAE_data', ...
    'ImageCorrelation_data', ...
    'EdgeDice_data', ...
    'GradientEnergy_data', ...
    'GradientRetention_data', ...
    'LaplacianVariance_data', ...
    'LaplacianRetention_data', ...
    'Entropy_data', ...
    'EntropyRetention_data', ...
    'HFPower_data', ...
    'HFRetention_data', ...
    'GLCMContrast_data', ...
    'GLCMCorrelation_data', ...
    'GLCMEnergy_data', ...
    'GLCMHomogeneity_data', ...
    '-v7.3');


fprintf('\n');
fprintf('============================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('Total images = %d\n',nFiles);
fprintf('Metric matrix size = %d x %d\n', ...
    nFiles,nResolution);
fprintf('Saved to:\n%s\n',outputFile);
fprintf('============================================\n');

end


%% =========================================================
% FUNCTION 1
% MAKE RGB
%% =========================================================

function I = makeRGB(I)

if ndims(I) == 2

    I = repmat(I,[1 1 3]);

elseif size(I,3) > 3

    I = I(:,:,1:3);

end

end


%% =========================================================
% FUNCTION 2
% FIND CORRESPONDING MASK
%% =========================================================

function maskFile = findMaskFile(maskFolder,imageName)

maskFile = '';

[~,baseName,~] = fileparts(imageName);


extensions = { ...
    '.png', ...
    '.tif', ...
    '.tiff', ...
    '.jpg', ...
    '.jpeg'};


possibleNames = { ...
    baseName, ...
    [baseName '_mask'], ...
    [baseName '-mask'], ...
    [baseName '_annotation']};


for a = 1:length(possibleNames)

    for b = 1:length(extensions)

        candidate = ...
            fullfile( ...
            maskFolder, ...
            [possibleNames{a} extensions{b}]);


        if exist(candidate,'file')

            maskFile = candidate;

            return;

        end

    end

end

end


%% =========================================================
% FUNCTION 3
% MASKED SSIM
%% =========================================================

function value = calculateMaskedSSIM(A,B,mask)

x = double(A(mask));

y = double(B(mask));


if isempty(x)

    value = NaN;

    return;

end


muX = mean(x);

muY = mean(y);


xZero = x - muX;

yZero = y - muY;


varianceX = ...
    mean(xZero.^2);


varianceY = ...
    mean(yZero.^2);


covarianceXY = ...
    mean(xZero .* yZero);


C1 = 0.01^2;

C2 = 0.03^2;


numerator = ...
    (2*muX*muY + C1) * ...
    (2*covarianceXY + C2);


denominator = ...
    (muX^2 + muY^2 + C1) * ...
    (varianceX + varianceY + C2);


if denominator ~= 0

    value = ...
        numerator / denominator;

else

    value = NaN;

end

end


%% =========================================================
% FUNCTION 4
% MASKED ENTROPY
%% =========================================================

function H = calculateMaskedEntropy(I,mask)

values = ...
    double(I(mask));


if isempty(values)

    H = NaN;

    return;

end


values(values < 0) = 0;

values(values > 1) = 1;


indices = ...
    floor(values*255) + 1;


counts = ...
    accumarray( ...
    indices(:), ...
    1, ...
    [256 1]);


totalCount = ...
    sum(counts);


if totalCount == 0

    H = NaN;

    return;

end


probability = ...
    counts / totalCount;


probability = ...
    probability(probability > 0);


H = ...
    -sum( ...
    probability .* ...
    log2(probability));

end


%% =========================================================
% FUNCTION 5
% MASKED HIGH-FREQUENCY POWER
%% =========================================================

function fraction = calculateMaskedHFPower(I,mask)

I = double(I);


[rowIndex,columnIndex] = ...
    find(mask);


if isempty(rowIndex)

    fraction = NaN;

    return;

end


row1 = min(rowIndex);

row2 = max(rowIndex);

column1 = min(columnIndex);

column2 = max(columnIndex);


Icrop = ...
    I(row1:row2,column1:column2);


Mcrop = ...
    mask(row1:row2,column1:column2);


tissueValues = ...
    Icrop(Mcrop);


if isempty(tissueValues)

    fraction = NaN;

    return;

end


%% Replace background by mean tissue intensity

meanTissue = ...
    mean(tissueValues);


Icrop(~Mcrop) = ...
    meanTissue;


%% Remove mean intensity

Icrop = ...
    Icrop - mean(Icrop(:));


%% Fourier transform

F = ...
    fftshift(fft2(Icrop));


P = ...
    abs(F).^2;


[rows,columns] = ...
    size(Icrop);


x = ...
    (1:columns) - ...
    (columns+1)/2;


y = ...
    (1:rows) - ...
    (rows+1)/2;


[X,Y] = ...
    meshgrid(x,y);


radius = ...
    sqrt(X.^2 + Y.^2);


maximumRadius = ...
    min(rows,columns)/2;


%% Upper 50% of spatial frequencies

highMask = ...
    radius >= ...
    (0.5*maximumRadius);


totalPower = ...
    sum(P(:));


highValues = ...
    P(highMask);


highPower = ...
    sum(highValues(:));


if totalPower > 0

    fraction = ...
        highPower / totalPower;

else

    fraction = 0;

end

end


%% =========================================================
% FUNCTION 6
% MASKED GLCM
%% =========================================================

function features = ...
    calculateMaskedGLCM(I,mask,numLevels)


I = double(I);


I(I < 0) = 0;

I(I > 1) = 1;


%% Quantize

Q = ...
    floor(I*numLevels) + 1;


Q(Q < 1) = 1;

Q(Q > numLevels) = numLevels;


%% Four directions

offsets = [ ...
     0  1;
    -1  1;
    -1  0;
    -1 -1];


[rows,columns] = ...
    size(I);


G = ...
    zeros(numLevels,numLevels);


for k = 1:size(offsets,1)

    dr = ...
        offsets(k,1);

    dc = ...
        offsets(k,2);


    rowStart = ...
        max(1,1-dr);

    rowEnd = ...
        min(rows,rows-dr);


    colStart = ...
        max(1,1-dc);

    colEnd = ...
        min(columns,columns-dc);


    rowA = ...
        rowStart:rowEnd;

    colA = ...
        colStart:colEnd;


    rowB = ...
        rowA + dr;

    colB = ...
        colA + dc;


    A = ...
        Q(rowA,colA);

    B = ...
        Q(rowB,colB);


    maskA = ...
        mask(rowA,colA);

    maskB = ...
        mask(rowB,colB);


    validPair = ...
        maskA & maskB;


    valuesA = ...
        A(validPair);

    valuesB = ...
        B(validPair);


    if isempty(valuesA)

        continue;

    end


    linearIndex = ...
        sub2ind( ...
        [numLevels numLevels], ...
        valuesA(:), ...
        valuesB(:));


    counts = ...
        accumarray( ...
        linearIndex, ...
        1, ...
        [numLevels*numLevels 1]);


    Gtemp = ...
        reshape( ...
        counts, ...
        numLevels, ...
        numLevels);


    G = ...
        G + Gtemp + Gtemp';

end


total = ...
    sum(G(:));


if total == 0

    features.Contrast = NaN;

    features.Correlation = NaN;

    features.Energy = NaN;

    features.Homogeneity = NaN;

    return;

end


P = ...
    G / total;


[ii,jj] = ...
    ndgrid( ...
    1:numLevels, ...
    1:numLevels);


%% Contrast

temp = ...
    ((ii-jj).^2) .* P;


features.Contrast = ...
    sum(temp(:));


%% Energy

temp = ...
    P.^2;


features.Energy = ...
    sum(temp(:));


%% Homogeneity

temp = ...
    P ./ ...
    (1 + abs(ii-jj));


features.Homogeneity = ...
    sum(temp(:));


%% Correlation

tempI = ...
    ii .* P;


tempJ = ...
    jj .* P;


muI = ...
    sum(tempI(:));


muJ = ...
    sum(tempJ(:));


temp = ...
    ((ii-muI).^2) .* P;


sigmaI = ...
    sqrt(sum(temp(:)));


temp = ...
    ((jj-muJ).^2) .* P;


sigmaJ = ...
    sqrt(sum(temp(:)));


if sigmaI > 0 && sigmaJ > 0

    temp = ...
        (ii-muI) .* ...
        (jj-muJ) .* ...
        P;


    numerator = ...
        sum(temp(:));


    features.Correlation = ...
        numerator / ...
        (sigmaI*sigmaJ);

else

    features.Correlation = NaN;

end

end