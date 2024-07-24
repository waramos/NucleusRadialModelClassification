function [PC, Clusterlabels] = CropPointCloud(PC, Clusterlabels, cropregion)
% CROPPOINTCLOUD will "crop" a point cloud to a region that the user 
% specifies as a 2x3 matrix where the first row is the starting index, the 
% second row is the ending index, and the columns correspond to dimensions
% x, y, z, respectively. This essentially throws out any points that are
% outside the requested region of interest.
%
% 
% INPUTS:
% 
% PC (Mx3 matrix) - Data struct that holds the original points from all 2D
% detections across the 3D volume.
%
% Clusterlabels (Mx1 vector) - the labels for the points in the
% point cloud that suggest what cluster points belong to. This will have
% any out of bound point dropped after bound correction
% 
% cropregion (2x3 matrix) - columns corresponding to lower and upper bound
% of the x, y, and z dimensions, respectively.
%
%
%
% OUTPUTS:
%
% PC (Mx3 matrix) - point cloud with ROI crop applied, i.e. points removed
% from set.
%
% Clusterlabels (Mx1 integer valued vector) - the labels for the points 
% in the point cloud that suggest what cluster points belong to. Will have
% any out of bound point dropped after bound correction
%
%
%
% William A. Ramos, Kumar Lab @ MBL July 2024

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

    % Adjust the point clouds by the shift
    PC = PC - [x1 y1 z1];
end