function [Fpatch,Hpatch] = ...
    extractPatchFeatures( ...
    FTIR,HE,tissueMask,...
    patchSize,minTissueFraction)

%% =========================================================
% FTIR
% 256 x 256 x 965
%
% H&E
% 256 x 256 x 3
%
% patch = 16 x 16
%% =========================================================

Fpatch = [];

Hpatch = [];


for r = 1:patchSize:256

    for c = 1:patchSize:256

        r2 = r + patchSize - 1;

        c2 = c + patchSize - 1;


        %% -------------------------------------------------
        % MASK FOR CURRENT PATCH
        %% -------------------------------------------------

        patchMask = ...
            tissueMask(r:r2,c:c2);


        tissueFraction = ...
            mean(patchMask(:));


        %% -------------------------------------------------
        % Reject background patches
        %% -------------------------------------------------

        if tissueFraction < minTissueFraction

            continue

        end


        %% =================================================
        % FTIR PATCH
        %% =================================================

        ftirPatch = ...
            FTIR(r:r2,c:c2,:);


        %
        % 16 x 16 x 965
        %     ->
        % 256 x 965
        %

        ftirPatch = ...
            reshape(ftirPatch,[],965);


        %
        % Only tissue pixels
        %

        ftirPatch = ...
            ftirPatch(patchMask(:),:);


        %
        % Mean spectrum
        %
        % 256 x 965
        %     ->
        % 1 x 965
        %

        meanSpectrum = ...
            mean(ftirPatch,1);


        %% =================================================
        % H&E PATCH
        %% =================================================

        hePatch = ...
            HE(r:r2,c:c2,:);


        %
        % 16 x 16 x 3
        %       ->
        % 1 x 768
        %

        heVector = ...
            reshape(hePatch,1,[]);


        %% =================================================
        % STORE MATCHED PATCH
        %% =================================================

        Fpatch = ...
            [Fpatch; meanSpectrum];

        Hpatch = ...
            [Hpatch; heVector];

    end

end

end