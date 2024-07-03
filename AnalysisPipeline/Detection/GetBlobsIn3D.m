function [C, Mask, P] = GetBlobsIn3D(I, s1, s2, thresh, r)
% FINDBLOBS will return both the reduced cluster points and image mask
% based on a difference of gaussian blob detection approach.

    % DoG
    C  = [];
    nz = size(I, 3);
    for i = 1:nz
        Im   = I(:,:,i);
        Mask = SimpleDoG(Im, s1, s2);
    
        % Threshold on DoG Image with morphological filtering and flood fill
        Mask = ThresholdCleanly(Mask, thresh, r);


        % Centroids of original detection
        Cn = regionprops(Mask, 'Centroid');
        Cn = cat(1, Cn.Centroid);
    
        % Get watershedded points
        P = SkeletonBasedWaterShed(Mask);

        % Finding centroids' nearest neighbor
        [~, D] = knnsearch(Cn, Cn, 'K', 2);
        C_l2   = D(:,2);
        C_l2   = floor(C_l2);
    
        % Radius approximation is average l2 norm
        r    = mean(C_l2)/2;
        rmsg = ['Mean radius: ' num2str(r)];
        disp(rmsg)
    
        % Drop points not within r of a centroid
        [~, D]   = knnsearch(Cn, P, 'K', 1);
        idx      = D<=r;
        [row, ~] = find(idx);
        P        = P(row, :);

        % Update on dropped points
        npdp = sum(~idx);
        dmsg = ['Num nuc. dropped: ' num2str(npdp)];
        disp(dmsg)

        % P = ClusterNucleiCentroids(P, r, CID, voxelsz);
        z = repmat(i, [size(P,1) 1]);
        P = [P z];
        C = cat(1, C, P);
        pd = (i/nz)*100;
        pd = [num2str(pd) ' % done'];
        disp(pd)
    end

    disp('Computing DBSCAN in 3D')
    CID     = 0.8;
    voxelsz = 1;
    P       = ClusterNucleiCentroids(C, r/2, CID, voxelsz);

end