clc
clear
close all

%% =========================================================
% INPUT FOLDER
%% =========================================================

folderPath = ...
'D:\Prostate Cancer Assignment\aligned_he_annotation_only';


%% =========================================================
% OUTPUT ROOT FOLDER
%% =========================================================

outputRoot = ...
'D:\Prostate Cancer Assignment\HE_spatial_resolution';

if ~exist(outputRoot,'dir')
    mkdir(outputRoot);
end


%% =========================================================
% Spatial resolutions
%% =========================================================

resolutionList = [1024 512 256 128 64];


%% =========================================================
% Find image files
%% =========================================================

files = [ ...
    dir(fullfile(folderPath,'*.png')); ...
  
    ];

fprintf('Total H&E images found = %d\n',numel(files));


%% =========================================================
% CREATE OUTPUT FOLDERS
%
% Native resolution folders:
%
% Resolution_1024
% Resolution_512
% Resolution_256
% Resolution_128
% Resolution_64
% Resolution_32
%
% Also create comparison folders where all reduced
% images are resized back to 1024 x 1024.
%% =========================================================

for r = resolutionList

    nativeFolder = fullfile( ...
        outputRoot, ...
        ['Resolution_' num2str(r)]);

    if ~exist(nativeFolder,'dir')
        mkdir(nativeFolder);
    end

end


%% Comparison folders

comparisonRoot = fullfile( ...
    outputRoot, ...
    'Upsampled_to_1024');

if ~exist(comparisonRoot,'dir')
    mkdir(comparisonRoot);
end


reducedResolutions = [512 256 128 64];

for r = reducedResolutions

    comparisonFolder = fullfile( ...
        comparisonRoot, ...
        ['From_' num2str(r)]);

    if ~exist(comparisonFolder,'dir')
        mkdir(comparisonFolder);
    end

end


%% =========================================================
% PROCESS ALL IMAGES
%% =========================================================

for i = 1:numel(files)

    fprintf('\n========================================\n');
    fprintf('Processing %d / %d\n',i,numel(files));
    fprintf('%s\n',files(i).name);
    fprintf('========================================\n');


    %% -----------------------------------------------------
    % Load original H&E image
    %% -----------------------------------------------------

    inputFile = fullfile( ...
        folderPath, ...
        files(i).name);

    I = imread(inputFile);


    %% -----------------------------------------------------
    % Ensure RGB
    %% -----------------------------------------------------

    if ndims(I) == 2

        I = repmat(I,[1 1 3]);

    elseif size(I,3) > 3

        I = I(:,:,1:3);

    end


    %% -----------------------------------------------------
    % Check original size
    %% -----------------------------------------------------

    fprintf('Original size = %d x %d x %d\n', ...
        size(I,1),size(I,2),size(I,3));


    %% -----------------------------------------------------
    % If original image is not exactly 1024 x 1024,
    % resize to reference size.
    %
    % If all images are already 1024 x 1024, this does
    % nothing.
    %% -----------------------------------------------------

    if size(I,1) ~= 1024 || size(I,2) ~= 1024

        warning(['Image %s is not 1024x1024. ' ...
                 'Resizing to reference size.'], ...
                 files(i).name);

        I = imresize( ...
            I, ...
            [1024 1024], ...
            'bicubic', ...
            'Antialiasing',true);

    end


    %% =====================================================
    % SAVE FULL-RESOLUTION REFERENCE
    %
    % Keep exactly the same filename.
    %% =====================================================

    referenceFile = fullfile( ...
        outputRoot, ...
        'Resolution_1024', ...
        files(i).name);

    imwrite(I,referenceFile);


    %% =====================================================
    % CREATE REDUCED SPATIAL RESOLUTIONS
    %% =====================================================

    for j = 1:length(reducedResolutions)

        R = reducedResolutions(j);


        fprintf('   Creating %d x %d\n',R,R);


        %% -------------------------------------------------
        % Downsample with anti-aliasing
        %
        % Field of view stays identical.
        %% -------------------------------------------------

        I_reduced = imresize( ...
            I, ...
            [R R], ...
            'bicubic', ...
            'Antialiasing',true);


        %% -------------------------------------------------
        % Save native reduced-resolution image
        %
        % KEEP SAME ORIGINAL FILENAME
        %% -------------------------------------------------

        nativeFolder = fullfile( ...
            outputRoot, ...
            ['Resolution_' num2str(R)]);


        nativeFile = fullfile( ...
            nativeFolder, ...
            files(i).name);


        imwrite(I_reduced,nativeFile);


        %% =================================================
        % RESIZE BACK TO 1024 x 1024
        %
        % This image should ONLY be used for comparison
        % against the original using:
        %
        % SSIM
        % PSNR
        % RMSE
        % edge preservation
        %
        % Upsampling DOES NOT restore lost information.
        %% =================================================

        I_compare = imresize( ...
            I_reduced, ...
            [1024 1024], ...
            'bicubic', ...
            'Antialiasing',true);


        %% -------------------------------------------------
        % Save comparison image
        %
        % Same filename because it is inside its own folder.
        %% -------------------------------------------------

        comparisonFolder = fullfile( ...
            comparisonRoot, ...
            ['From_' num2str(R)]);


        comparisonFile = fullfile( ...
            comparisonFolder, ...
            files(i).name);


        imwrite(I_compare,comparisonFile);

    end


    fprintf('Finished: %s\n',files(i).name);

end


fprintf('\n========================================\n');
fprintf('ALL H&E IMAGES COMPLETED\n');
fprintf('========================================\n');