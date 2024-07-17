function [P, refidx] = MaskContainsMask(Refmask, Mask)
% MASKCONTAINSMASK will return a T/F classification regarding whether a
% reference mask contains some overlap with another mask.

    % Voxel indices for masks
    refidx  = regionprops3(Refmask, 'VoxelIdxList'); % refernce mask

    % Intersection of the two
    I       = Refmask | Mask;
    intidx  = regionprops3(I, 'VoxelIdxList');  % voxel indices of intersecting regions
    allint  = [intidx.VoxelIdxList];

    % Finding the original masked ROIs with the overlapping regions
    numROIs = size(refidx, 1);
    P       = false(numROIs, 1);
    for i = 1:numROIs
        % Pulling out reference indices for a given ROI
        idx  = refidx(i).VoxelIdxList;
        o    = intersect(idx, allint);
        P(i) = isempty(o);
    end
end