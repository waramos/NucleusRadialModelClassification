%% Analysis of Overlap with Radial Model
% This code aims to compare a 3D pointcloud from detected nuclei to a stack
% of segmented slices from an image volume to analyze any overlap that
% might be present in the segmented channel. Nuclei are then classified as
% either positive or negative for overlap with the segmented signal. The 
% former should be double precision as a 3-column array with N-rows for N 
% number of points. The stack of segmented slices is to be an M x N x Z 
% logical array with M rows, N columns, and Z slices.


%% Variables to consider
% P is the pointcloud representing the nuclei detections in 3D
% r is the radius of the spherical model in pixels/voxels
% V is the binary volume (logical array, i.e., segmented stack)


%% Simulated Results
% Considering the anisotropic resolution, it might be good to interpolate
% the masks and then scale the points so that they are on the same
% coordinate system.
% e.g. 
% if the camera pixel size is 6.5 µm, objective magnifaction is 100 and x/y
% pixels have the same resolution but the image volume was taken with 
% 0.3 µm z spacing between planes:
% dx = dy = 6.5/100
% dz = 0.3, so then our z slices are spaced 4.6 pixels apart
% This can be taken into consideration later

%% Parameters
% Anisotropic Resolution
% Voxel spacing - 65 nm xy spacing, 300 nm z spacing
dx = 6.5/100; 
dy = dx;
dz = 0.3;

% Radius of nuclear radial model: 
% Assuming 3~4 µm diameter, we might want to give a micron beyond the
% nucleus size: e.g. 5-6 would be reasonable diameters
d = 6;
r = d/(2*dx); % Converts diameter to radius in pixels


%% Data Import
% Place holder for: Pointcloud - Detect Nuclei Centroids
P_x = randi(M*100, [100 1])/100;                            % CHANGE THIS
P_y = randi(N*100, [100 1])/100;
P_z = randi(Z*100, [100 1])/100;
P   = [P_x' P_y' P_z']; % Point cloud

% Place holder for: Binary Stack - Segmented Overlap Channel
Maskstack = zeros(2048, 3848, 20);                          % CHANGE THIS


%% Computation
% Finding Nuclei with overlap
[Nuc_pos,...
 Nuc_all,...
 CCData,...
 localizationerror] = RadialOverlap(P, r, Maskstack, [dx dy dz]);


%% Visualization
[x_, y_, z_] = ind2sub(size(Maskstack), idx);
AS           = alphaShape(x_, y_, z_); 
hold on
plot(AS, 'FaceColor', 'magenta', 'FaceAlpha', 0.15, 'LineStyle', 'none')