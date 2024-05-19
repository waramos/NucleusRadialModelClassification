function Voxcount = ROIInSpheres(C, Mask, r)
% ROIINSPHERES will count the number of voxels from the 3D mask volume that
% are within radius, r, of points in the pointcloud P. Beyond a simple
% classification, this will enable some statistical measurement as well.

    % Cleanup mask and find pixels
    Mask  = imerode(Mask, [0 1 0; 1 1 1; 0 1 0]);
    idx   = find(Mask);
    [X,...
     Y,...
     Z]   = ind2sub(size(Mask), idx);
    P     = [X Y Z];

    % Update on voxels found
    prct = sum(Mask(:)==1)/numel(Mask);
    msg  = ['Voxel indices found. ' newline ...
            num2str(prct*100, '%.2f') '% voxels inside ROI'];
    disp(msg)

    % Checks number of points and initializes a voxel count vector
    numpts   = size(C, 1);
    Voxcount = zeros(numpts, 1);
    for i = 1:numpts
        % Find number of voxels within radius, r, of centroid
        pt          = C(i,:);           % Point to compare against
        D           = pdist2(pt, P);
        D           = D<=r;
        Voxcount(i) = sum(D);  % Number of voxels within distance r
        prct        = (i/numpts)*100;
        msg         = ['Centroids analyzed: ' num2str(prct, '%.2f') '%'];
        disp(msg)
    end
end