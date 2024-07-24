


rsfactor   = dx/dz;
voxthresh  = 49;
visnhoods  = false;
savenhoods = false;
Results    = [];

for ch = 1:nchs
    for f = 1:nfiles
        % Segmentation of the channel of interest 
        % Pulling image to analyze
        Im     = I(:, :, :, f, ch);
        ChMask = HardThreshold(Im, thresh); % Segmentation
        ChMask = RefineMask(ChMask, 3);     % Cleans up segmentation

        if ~all((diff(cropregion,1) + 1) == ([sz(2) sz(1) sz(3)]))
            ChMask     = CropVol(ChMask, cropregion);
        end
        
        % Size information
        sz = size(ChMask, [1 2 3]);

        % Performs the classification and measurement
        [Nuclei_mask,...
         PC,...
         Clusters] = NuclearPointsClouds2Masks(Nuclei, f, r, sz, cropregion, rsfactor);
        NewResults = ClassifyExpressionFast(PC,Clusters, Nuclei_mask, ChMask, I(:,:,:,f, ch), r, voxthresh, true, savenhoods);
        Results = cat(1, Results, NewResults);
    end
    AllResults{ch} = Results;
end

