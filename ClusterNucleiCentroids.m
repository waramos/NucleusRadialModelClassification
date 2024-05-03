function [centroids, corepts, negpts, discardpercentage] = ClusterNucleiCentroids(P, r, CID, voxelsz)
% CLUSTERNUCLEICENTROIDS will use DBSCAN with a given radius, r to cluster
% detections from an image volume. It is assumed detections were done in 2D
% and the input is a pointcloud, P. Voxel spacing is used for determining
% minimum number of points required for a valid cluster/nucleus, given that
% the user has a minimum confidence requirement.

    % unpacking voxel dims
    dx = voxelsz(1);
    dy = voxelsz(2);
    dz = voxelsz(3);

    % Min num. points / slices for a detection, given confidence thresh
    rpx    = r/dz;            % radius in pixels
    minpts = ceil(2*rpx*CID); % minimum number of points

    % DBSCAN applied to points
    [idx, corepts] = dbscan(P, r, minpts);

    % If negative points and discard % requested
    if nargout >= 3
        negidx            = idx == -1;
        [neg_ridx, ~]     = find(negidx);
        negpts            = P(neg_ridx, :);
        discardpercentage = sum(negidx)/numel(idx);
    end
    
    % Initializing centroids
    cidx        = idx(idx~=-1);        % only grab positive detections
    centroidIdx = unique(cidx);        % unique centroids
    numC        = numel(centroidIdx);  % number of centroids
    centroids   = zeros(numC, 3);      % centroids array

    for i = 1:numC
        % Centroid index value
        ind            = cidx(i);
        [rowidx, ~]    = find(idx == ind);
        pt             = P(rowidx, :);
        centroids(i,:) = mean(pt, 1);
    end
end