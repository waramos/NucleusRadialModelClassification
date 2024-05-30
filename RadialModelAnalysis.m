%% Radial Distance Model: Analysis of Microscopy Image Volumes
% This code aims to compare a 3D pointcloud from detected nuclei to a stack
% of segmented slices from an image volume to analyze any overlap that
% might be present in the segmented channel. Nuclei are then classified as
% either positive or negative

% Example for the analysis:
% First rough pass of Michael Palmer's data with nuclei and two different
% channels.

OverallTimer = tic;

%% Loading in Nuclei pointclouds
% This is data exported from the Segmentation GUI as a point cloud from
% blob and/or puncta detections

% Location of data
datafolder = 
Nuc_file   = 'well1_2024_Feb_29_15.40.31.112_405_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat';
fid        = fullfile(datafolder, Nuc_file);

% Only read in the table with coordinates
load(fid, 'ResultsTable')

% Imaging parameters
lensmag      = 25;              % Lens magnification
campxsz      = 6.5;             % camera pixel size
xypxsz       = campxsz/lensmag; % xy pixel size
zspacing     = 0.3;             % z spacing in microns
ispixelunits = true;            % point cloud left in pixels for now

% Converting datatable to an isotropic point cloud of nuclei detections
P = SetPointsCoordinateSys(ResultsTable, xypxsz, zspacing, ispixelunits);


% CL update 1
msg = 'Step 1: Point cloud loaded succesfully';
disp(msg)

%% Nuclei Centroids from Clustering
% DBSCAN used here to consider tossing out "noisy" detections. We can be
% confident that if a nucleus appears on multiple slices, you should expect
% 1 detection per nucleus on each slice. A nucleus we can be confident in
% counting should be more than half inside the image volume, given that the
% centroid will be skewed away from its actual true location in a detection
% if too few of slices are used to compute it.

% Nucleus params
nuc_d   = 12;      % nucleus diameter
r       = nuc_d/2; % nucleus radius in microns
CID     = 0.75;     % Confidence In Detection - 75%

% Whether to use pixels or microns
if ispixelunits
    rspacing = 1;
else
    rspacing = xypxsz;
end

% Cluster nuclei and finds their centroids w/ detection confidence thresh
[centroids, corepts, negpts, discardpercentage] = ClusterNucleiCentroids(P, r, CID, rspacing);


% CL update 2
msg = 'Step 2: Point clouds clustered into centroids succesfully';
disp(msg)

%% Loading in Segmented Channels
fid     = fullfile(datafolder, "well1_2024_Feb_29_15.39.57.846_488_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat");
V1      = ExtractVolumeResults(fid);
Ch1Name = 'Col2a1';

fid     = fullfile(datafolder, "well1_2024_Feb_29_15.39.24.566_640_0.20502ms_2_Segmentation_FirstRoughPassInfo.mat");
V2      = ExtractVolumeResults(fid);
Ch2Name = 'ColX';


% CL update 3
msg = 'Step 3: Segmented channels loaded succesfully';
disp(msg)

%% Resampling Segmented Masks
% Resampling will take place to ensure that the segmented masks are also
% isotropically spaced.

% Timer
t = tic;

if ispixelunits
    dx = 1;
    dz = zspacing/xypxsz;
else
    dx = xypxsz;
    dz = zspacing;
end

% Isotropic Resample
V1 = IsoResampleVolume(V1, dx, dz);
V2 = IsoResampleVolume(V2, dx, dz);

% CL update 4
msg = 'Step 4: Binary volumes resampled succesfully';
disp(msg)

% Timer
T   = toc(t);
msg = ['Time to resample binary volumes:' newline ...
       num2str(T) ' seconds.'];
disp(msg)

%% Voxels inside Radius
% Finds number of voxels that are within the radius, r, of the nuclear
% centroids. Classification can then be done with some user selected
% threshold for voxel count; considering that there is a variety of sizes
% (i.e. number of voxels) across the 

% Timer
t = tic;

Voxelscount1 = ROIInSpheres(centroids, V1, r);
Voxelscount2 = ROIInSpheres(centroids, V2, r);

% CL update 5
msg = 'Step 5: Voxels within distance found succesfully';
disp(msg)

% Timer
T   = toc(t);
msg = ['Time to find voxels within distance r:' newline ...
       num2str(T) ' seconds.'];
disp(msg)

%% Covariance and Dependency Analysis
% We can consider how the two channels are related by analyzing their
% covariance as well as how they soread in relation to each other by
% computing the eigen vectors of the two channels by plotting them 
figure


CV = cov(Voxelscount1, Voxelscount2);

STD1 = std(Voxelscount1(:));
STD2 = std(Voxelscount2(:));
CorrCoeff = CV./[STD1^2 STD1*STD2; STD1*STD2 STD2^2];

hHM = heatmap(CV, ...
      "Title", 'Voxel Count Covariance', ...
      'XLabel', 'Variance of Channel', ...
      'YLabel', 'Variance of Channel');

% Covariance Labels
cv1name = ['\sigma(' Ch1Name ')'];
cv2name = ['\sigma(' Ch2Name ')'];
hHM.XDisplayLabels = {cv1name, cv2name};
hHM.YDisplayLabels = {cv1name, cv2name};


% PCA analysis
[Mj, Mn] = MomentsOfInertia(Voxelscount1, Voxelscount2);

BothChCounts = [Voxelscount1, Voxelscount2];
[pca_coeff, pca_score, pca_latent] = pca(BothChCounts);

% Plotting the moments of inertia for the data (e.g. spread)
figure
plot(Mj(:,1), Mj(:,2))
hold on

plot(Mn(:,1), Mn(:,2))

% 0 lines
xline(0)
yline(0)

plot(Voxelscount1, Voxelscount2, '.r', 'MarkerSize', 12);
axis equal


xlabel([Ch1Name ' Voxel Counts']); 
ylabel([Ch2Name ' Voxel Counts']);
ttext = ['Voxels Within ' num2str(r) 'µm of Nucleus Centroid' newline ...
         '(Data Spread with Moments of Inertia)'];
title(ttext)

% CL update 6
msg = 'Step 6: Covariance and PCA analysis complete';
disp(msg)

%% Thresholding Voxel Counts
% Here, we find all connected components with connectivity of 26.
% We then find the smallest connected component and assume this is the
% minimum number of voxels that would need to be found within the radial
% distance described above to count as a 

% Grab CC's and find smallest object
Ch1CC        = bwconncomp(V1);
Ch1CCVoxList = [Ch1CC.PixelIdxList];
Ch1CCVoxList = cellfun(@numel, Ch1CCVoxList);
Ch1VoxThresh = min(Ch1CCVoxList);

Ch2CC        = bwconncomp(V2);
Ch2CCVoxList = [Ch2CC.PixelIdxList];
Ch2CCVoxList = cellfun(@numel, Ch2CCVoxList);
Ch2VoxThresh = min(Ch2CCVoxList);

% CL update 7
msg = 'Step 7: Found smallest CCs per channel';
disp(msg)

%% Nucleus Classification
% The simple classification applies the threshold for the total number of
% voxel counts in the smallest object that was segmented for a given
% channel. If a nucleus has more voxel counts than the threshold, it is
% considered positive for that channel.

% Number of nuclei
numNuclei = size(centroids, 1);

% Positive classification per channel
Ch1Class = Voxelscount1 > Ch1VoxThresh;
Ch2Class = Voxelscount2 > Ch2VoxThresh;

% Consider Nuclei that have both or none
Ch1Only   = Ch1Class & ~BothChPos;
Ch2Only   = Ch2Class & ~BothChPos;
BothChPos = Ch1Class & Ch2Class;
NoCh      = ~Ch1Class & ~Ch2Class;

% Total number of nuclei that contain sufficient counts
Ch1Total  = sum(Ch1Only);
Ch2Total  = sum(Ch2Only);
BothTotal = sum(BothChPos);
NoneTotal = sum(NoCh);
Totals    = [NoneTotal; Ch1Total; Ch2Total; BothTotal];
Totals    = Totals/numNuclei;

% Alternatively, we could show this as a radial distriubution and draw on
% inspiration from statistical mechanics.. The interpretation then would
% be, what is the probability that you find a given collagen type or
% protein a certain distance away from the nucleus [centroid/center].
% However, this computation... at least at a glance, would be
% computationally expensive.


% CL update 8
msg = 'Step 8: Nuclei classified succesfully';
disp(msg)

%% Simple Stacked Bar Plots
% Will show proportion of nuclei that fall in each category and then how
% positive nuclei fall in different 

figure
varnames = {'None', Ch1Name, Ch2Name, 'Both'};
bar(1, Totals, 'stacked')
title('Nuclei')
ylabel('% of Nuclei')
legend(varnames)
xticklabels('Control')


%% Positive Nuclei Volume Centroids

% Ch 1 positive nuclei
[r, ~]       = find(Ch1Class);
positiveNuc1 = centroids(r, :);

% Ch 2 positive nuclei
[r, ~]       = find(Ch2Class);
positiveNuc2 = centroids(r, :);

%% Visualization

% Plotting Timer
t = tic;

% Figure and axes
[M, N, Z] = size(V1);
f         = figure;
ax        = axes(f);

% colors / colormap
c0 = [1 0.1 0.1];
c1 = [0.3 0.3 1];

% Plotting ALL nuclei
p0 = plot3(centroids(:,1), centroids(:,2), centroids(:,3), 'Color', c0, 'Marker', '.', 'LineStyle', 'none');
hold on

axis equal

% Plotting positive nuclei for...
% Channel 1
x    = positiveNuc1(:,1);
y    = positiveNuc1(:,2);
z    = positiveNuc1(:,3);
p0_1 = plot3(x, y, z, 'Color', fliplr(c0), 'Marker', 'o', 'LineStyle', 'none');    % Ch1 positive

% Channel 2
x    = positiveNuc2(:,1);
y    = positiveNuc2(:,2);
z    = positiveNuc2(:,3);
p0_2 = plot3(x, y, z, 'Color', fliplr(1-c0), 'Marker', 'o', 'LineStyle', 'none');  % Ch2 positive


% Timer
T   = toc(t);
msg = ['Time to plot:' newline ...
       num2str(T) ' seconds.'];
disp(msg)

%% Alpha shape visualization
% This visualization then shows the surface of the segmented regions

% % Channel 1 Signal 3D Alpha Shape
% idx1            = find(V1);
% [x_1, y_1, z_1] = ind2sub([M N Z], idx1);
% AS1             = alphaShape(x_1, y_1, z_1);
% 
% % Channel 2 Signal 3D Alpha Shape
% idx2            = find(V2);
% [x_2, y_2, z_2] = ind2sub([M N Z], idx2);
% AS2             = alphaShape(x_2, y_2, z_2);
% 
% % Plotting the two channels
% p1 = plot(AS1, 'FaceColor', c1, 'FaceAlpha', 0.15, 'LineStyle', 'none');
% p2 = plot(AS1, 'FaceColor', fliplr(c1), 'FaceAlpha', 0.15, 'LineStyle', 'none');
% 
% % Auto setting to full FOV w/ isotropic scale
% axis(ax, 'equal')
% ax.YLim   = [1 M];
% ax.XLim   = [1 N];
% ax.ZLim   = [1 Z];


% Timer
T   = toc(OverallTimer);
msg = ['Overall load and compute time:' newline ...
       num2str(T) ' seconds.'];
disp(msg)