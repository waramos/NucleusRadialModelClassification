function [Nucleimask, PC, Clusters] = NuclearPointsClouds2Masks(Nucleistruct, fidx, r, sz, cropregion, rsfactor)
% NUCLEARPOINTCLOUDS2MASKS will grab the datastruct that is stored in a
% NuclearDetections.mat file. That file has a datastruct that results from
% the DBSCAN processing of data from our segmentation GUI that uses a
% well known blob detection algorithm, Difference of Gaussians, to segment
% cell nuclei.
% It is expected that the resulting data struct will at least have the
% OriginalPoints field and the Cluster field.
% The OriginalPoints fields holds the 2D detections as a 3D pointcloud
% with z information based on the z plane that the detections originated
% from. The cluster field then has the label for each point in the point
% cloud, thus specifying which cluster (i.e. cell) it belongs to.
% 
%
% INPUTS:
% 
% Nucleistruct (1xM or Mx1 datastruct) - Data struct that holds the DBSCAN
% processed blob detection points, original points, and other information.
%
% fidx (scalar integer) - file index between 1 and M that will specify which
% set of nuclei detections you would like to look at.
% 
% sz (1x3 row vector) - size vector holds the dimension sizing of the image
% that the points will map to. This is needed for bound checking the
% indices of pixels that the points will map to
%
%
%
% OUTPUTS:
%
% Nuclei (logical 3D array) - binary image volume represents the nuclear
% detections across each 2D slice then morphologically dilated by r pixels
%
% PC (Mx3 matrix) - point cloud with crop and isotropic transformation
% applied along z
%
% Clusters (Mx1 integer valued vector) - the labels for the points in the
% point cloud that suggest what cluster points belong to. This will have
% any out of bound point dropped after bound correction
%
%
%
% William A. Ramos, Kumar Lab @ MBL July 2024


    % Default rescale factor assumes isotropic resolution
    if nargin < 6 || isempty(rsfactor)
        rsfactor = 1;
    end

    % Default "crop" grabs everything
    if nargin < 5 || isempty(cropregion)
        cropregion = [ones(1,3); sz];
    end
    
    % Pulling out original PC and the cluster labels
    PC       = Nucleistruct(fidx).OriginalPoints;
    Clusters = Nucleistruct(fidx).Cluster;
    
    % Corrections and cropping as needed - puts PC in image coords
    PC             = AdjustPointCloud(PC, sz, rsfactor);
    [PC, Clusters] = CropPointCloud(PC, Clusters, cropregion);
    
    % Binary image from nuclei detections
    Nucleimask = Points2Cylinders(PC, r, sz);

    % Only need to perform a crop when the size is not already correct,
    % otherwise, this implies the crop was done first for efficiency in the
    % case of batch processing data
    if ~all((diff(cropregion,1) + 1) == ([sz(2) sz(1) sz(3)]))
        Nucleimask = CropVol(Nucleimask, cropregion);
    end
end