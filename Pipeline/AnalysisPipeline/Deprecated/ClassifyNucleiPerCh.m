function Results = ClassifyNucleiPerCh(PC, Clusterlabels, r, Chmask, coordsys, voxelsize, cropregion)

    % Coordinate system
    if nargin < 4 || isempty(coordsys)
        coordsys = 'isotropic';
    end

    % Voxel size: [dx dy dz]
    if nargin < 5 || isempty(voxelsize)
        voxelsize = [1 1 1];
    end

    % Point clouds are to be scaled back to image coordinates
    switch coordsys
        case 'isotropic'
            dz_dx = voxelsize(3)/voxelsize(1);
        case 'image'
            dz_dx = 1;
    end

    % Isotropically rescales the third dimension of the point cloud
    PC(:,3) = PC(:,3)*dz_dx;

    % In case a crop was requested
    if ~isempty(cropregion)
        Chmask          = CropVol(Chmask, cropregion);
        [PC,...
         Clusterlabels] = CropPC(PC, Clusterlabels, cropregion);

        PC(:,1) = PC(:,1) - cropregion(1,1) + 1;
        PC(:,2) = PC(:,2) - cropregion(1,2) + 1;
        PC(:,3) = PC(:,3) - cropregion(1,3) + 1;
    end

    % Image bounds
    sz = size(Chmask, [1 2 3]);

    % Nuclei represented by cylinders
    Nuclei = Points2Cylinders(PC, r, sz);
    Results.Nuclei = Nuclei;

    % Converting binary volume to connected components w pixel list
    NucCC = bwconncomp(Nuclei, 8);
    NucRP = regionprops(NucCC);
    ChCC  = bwconncomp(Chmask);
    ChRP  = regionprops(ChCC);

    % Get overlapping mask
    O = Chmask & Nuclei;
    ORegions = bwconncomp(O, 8);
    OPixels = [ORegions.PixelIdxList];
    

end


function Vol = CropVol(Vol, cropregion)
% CROPVOL will crop a volume based on the region specified by the
% cropregion matrix. First row of matrix is starting index and second row
% is the ending index.

    % Get crop coordinates
    x   = cropregion(1,1):cropregion(2,1);
    y   = cropregion(1,2):cropregion(2,2);
    z   = cropregion(1,3):cropregion(2,3);
    Vol = Vol(y, x, z);
end


function [PC, Clusterlabels] = CropPC(PC, Clusterlabels, cropregion)
% CROPPC will "crop" a point cloud to a region that the user specifies as a
% 2x3 matrix where the first row is the starting index, the second row is
% the ending index, and the columns correspond to dimensions. This function
% also expects the user to pass in cluster labels so that the labels
% corresponding to remaining points in the point cloud can be preserved.

    % Getting bounds of crop region
    x1 = cropregion(1,1);
    x2 = cropregion(2,1);
    y1 = cropregion(1,2);
    y2 = cropregion(2,2);
    z1 = cropregion(1,3);
    z2 = cropregion(2,3);

    % Splitting dimensions
    X = PC(:,1);
    Y = PC(:,2);
    Z = PC(:,3);

    % Getting valid indices 
    idx = (X>x1) & (X<x2);
    idy = (Y>y1) & (Y<y2);
    idz = (Z>z1) & (Z<z2);

    % Getting valid rows
    ind  = idx & idy & idz;
    rows = find(ind);

    % Extract valid points and their labels
    PC            = PC(rows, :);
    Clusterlabels = Clusterlabels(rows, :);
end





function Nuclei = Points2Cylinders(PC, r, sz)

    % Round to nearest pixel - new var enables comparison
    PC_pix = round(PC);
    y      = PC_pix(:,1);
    x      = PC_pix(:,2);
    z      = PC_pix(:,3);

    % Convert point cloud to linear indices corresponding to pixels
    idx = sub2ind(sz, x, y, z);

    % Init nuclei array and set the true values
    Nuclei      = false(sz(1), sz(2), sz(3));
    Nuclei(idx) = true;

    % Convolve Nuclei stack with disk kernal for cylinder model
    h      = MakeDiskKernel(r);
    Nuclei = imdilate(Nuclei, h);
end


function h = MakeDiskKernel(r)
% MAKEDISKKERNEL will make a disk kernel based off of a given radius
    h = -r:r;
    h = h.^2;
    h = h+h';
    h = sqrt(h)<=r;
end


function R = RadialProbPerNucleus(pts, r, ch)
end