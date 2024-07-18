function P = RefinedBlobDetection(Mask)
% REFINEDBLOBDETECTION will find the centroids, C, of connected components
% (CCs) in a mask and then produce a watershedded version of the mask and
% get the centroids, P, of the resulting CCs. Any points from the set, P,
% further than distance, r, of any given centroid from C are discarded. The
% remaining points are considered valid detections to be let alone in each
% 2D slice that was segmented. This is done for every slice in the stack of
% masks. Note that distance, r, is automatically calculated as half of the
% mean distance between nearest neighbors.
% radius are included. The algorithm automatically determines the
% appropriate radius to use: it checks the nearest neighbor 
    C = regionprops(Mask, 'Centroid');
    C = cat(1, C.Centroid);

    % Clusters of points from watershed
    Mask = gpuArray(Mask);
    P    = SkeletonBasedWaterShed(Mask);

    % Finding centroids' nearest neighbor
    [~, D] = knnsearch(C, C, 'K', 2);
    C_l2   = D(:,2);
    C_l2   = floor(C_l2);

    % Radius approximation
    r = mean(C_l2)/2;

    % Drop points not within r of a centroid
    [~, D]   = knnsearch(C, P, 'K', 1);
    idx      = D<=r;
    [row, ~] = find(idx);
    P        = P(row, :);
end