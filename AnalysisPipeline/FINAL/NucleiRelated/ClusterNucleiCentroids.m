function [centroids, corepts, negpts, discardpercentage, idx] = ClusterNucleiCentroids(P, r, CID, voxelsz)
% CLUSTERNUCLEICENTROIDS will use DBSCAN with a given radius, r to cluster
% detections from an image volume. It is assumed detections were done in 2D
% and the input is a pointcloud, P. If r is not in pixel units, the voxel
% size can be given as well. 

    % Min num. points / slices for a detection, given confidence thresh
    if nargin<4 || ~isempty(voxelsz)
        rpx     = r;
    else
        rpx     = r/voxelsz(3);          % radius in pixels
        rfactor = voxelsz(3)/voxelsz(1); % in case of a lack of isotropic resolution
        rpx     = rpx * rfactor;
    end

    % Can force min number of points per cluster
    if nargin < 3 || isempty(CID)
        minpts = 1;
    else
        minpts = ceil(rpx*CID); % minimum number of points
    end

    % DBSCAN applied to points
    [idx, corepts] = dbscan(P, r, minpts);

    % If negative points and discard % requested
    if nargout >= 3
        negidx            = idx == -1;
        [neg_ridx, ~]     = find(negidx);
        negpts            = P(neg_ridx, :);
        discardpercentage = sum(negidx)/numel(idx);
    end
    
    % Initializing centroids array
    cidx        = idx(idx~=-1);        % only grab positive detections
    centroidIdx = unique(cidx);        % unique centroids
    numC        = numel(centroidIdx);  % number of centroids
    centroids   = zeros(numC, 3);      % centroids array

    for i = 1:numC
        % Centroid index value
        ind            = centroidIdx(i);
        [rowidx, ~]    = find(idx == ind);
        pt             = P(rowidx, :);
        centroids(i,:) = mean(pt, 1);
    end
end