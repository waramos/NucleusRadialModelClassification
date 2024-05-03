%% Radial Distance Model: Analysis of Microscopy Image Volumes
% This code aims to compare a 3D pointcloud from detected nuclei to a stack
% of segmented slices from an image volume to analyze any overlap that
% might be present in the segmented channel. Nuclei are then classified as
% either positive or negative

% Example for the analysis:
% First rough pass of Michael Palmer's data with nuclei and two different
% channels.

%% Loading in Nuclei pointclouds
% This is data exported from the Segmentation GUI as a point cloud from
% blob and/or puncta detections

% Location of data
datafolder = 'D:\2024_02_29 - Michael Palmer Cells\Analysis';
Nuc_file   = 'well1_2024_Feb_29_15.40.31.112_405_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat';
fid        = fullfile(datafolder, Nuc_file);

% Only read in the table with coordinates
load(fid, 'ResultsTable')

xypxsz       = 6.5/25; % xy pixel size
zspacing     = 0.3;    % z spacing in microns
ispixelunits = true;   % point cloud left in pixels for now

% Converting datatable to an isotropic point cloud of nuclei detections
P = SetPointsCoordinateSys(ResultsTable, xypxsz, zspacing, ispixelunits);

%% Nuclei Centroids from Clustering
% DBSCAN used here to consider tossing out "noisy" detections. We can be
% confident that if a nucleus appears on multiple slices, you should expect
% 1 detection per nucleus on each slice. A nucleus we can be confident in
% counting should be more than half inside the image volume, given that the
% centroid will be skewed away from its actual true location in a detection
% if too few of slices are used to compute it.

% Voxel sizing
dx      = 1;
dy      = 1;
dz      = zspacing/xypxsz;
voxelsz = [dx dy dz];

% Nucleus params
nuc_d   = 12;      % nucleus diameter
r       = nuc_d/2; % nucleus radius in microns
CID     = 0.9;     % Confidence In Detection - 90%

% Cluster nuclei and finds their centroids w/ detection confidence thresh
[centroids, corepts, negpts, discardpercentage] = ClusterNucleiCentroids(P, r, CID, voxelsz);

cpts = corepts(corepts==1);

%% error calculations

% Comparing number of cores vs number of centroids - should be the same
% numCentroids = numel(centroids);
% numCores     = numel(corepts);




%% Loading in Ch 1

fid = fullfile(datafolder, "well1_2024_Feb_29_15.39.57.846_488_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat");
V1  = ExtractVolumeResults(fid);

%% Loading in Ch 2

fid = fullfile(datafolder, "well1_2024_Feb_29_15.39.57.846_488_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat");
V2  = ExtractVolumeResults(fid);

%% Nuclei Classification
% Finding Nuclei with overlap

% Radius is in pixels here
[positiveNuc1,...
 allNuc1,...
 ccData1,...
 localizationMSE] = RadialOverlap(centroids, r/dz, V1, voxelsz, true);

% ALTERNATIVELY:
% Can likely do an inpolygon check
% 


%% Quantification and Stats
% Of the nuclei that we deemed to be sufficiently accurate detected, how
% many were classified as positive for channel 1?
FractionPos = (ccData1.Overlapping.NumObjects/ccData1.Original.NumObjects)*100;

%% Converting image back to points
[M, N, Zrs] = size(allNuc1, [1 2 3]);
idx_pos     = find(positiveNuc1(:));
[pos_x,...
 pos_y,...
 pos_z]     = ind2sub([M N Zrs], idx_pos); % Ch 1 positive nuclei

%% Visualization

% Figure and axes
[M, N, Z] = size(V1);
f         = figure;
ax        = axes(f);

% colors / colormap
c0 = [1 0.1 0.1];
c1 = [0.3 0.3 1];

% Plotting nuclei
p0 = plot3(P(:,1), P(:,2), P(:,3), 'Color', c0, 'Marker', '.', 'LineStyle', 'none');
hold on

% Consider radial model and neg/pos nuclei
% e.g.
p0_1 = plot3(P(:,1), P(:,2), P(:,3), 'Color', fliplr(c0), 'Marker', 'o', 'LineStyle', 'none');    % Ch1 positive
p0_2 = plot3(P(:,1), P(:,2), P(:,3), 'Color', fliplr(1-c0), 'Marker', 'o', 'LineStyle', 'none');  % Ch2 positive

% Channel 1
idx1            = find(V1);
[x_1, y_1, z_1] = ind2sub([M N Z], idx1);
AS1             = alphaShape(x_1, y_1, z_1);

% Channel 2
idx2            = find(V2);
[x_2, y_2, z_2] = ind2sub([M N Z], idx2);
AS2             = alphaShape(x_2, y_2, z_2);

% Plotting the two channels
p1 = plot(AS1, 'FaceColor', c1, 'FaceAlpha', 0.15, 'LineStyle', 'none');
p2 = plot(AS1, 'FaceColor', fliplr(c1), 'FaceAlpha', 0.15, 'LineStyle', 'none');

% Auto setting to full FOV w/ isotropic scale
axis(ax, 'equal')
ax.YLim   = [1 M];
ax.XLim   = [1 N];
ax.ZLim   = [1 Z];
