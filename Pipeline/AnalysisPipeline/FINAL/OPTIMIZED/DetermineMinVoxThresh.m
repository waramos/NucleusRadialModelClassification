function thresh = DetermineMinVoxThresh(Mask)
% DETERMINEMINVOXTHRESH determines the minimum number of voxels needed to
% be overlapping for a positive classification.
    CC         = bwconncomp(Mask, 8);
    RP         = regionprops(CC, 'Area');
    A          = [RP.Area];
    [N, edges] = histcounts(A, 1:1000);
    [~, idx]   = max(N);
    thresh     = edges(idx);
end