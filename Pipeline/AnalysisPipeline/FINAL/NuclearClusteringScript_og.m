
%% Loading in data
[fname, fpath] = uigetfile('*.mat', 'Select Nuclear Detection Data');
% addpath 'D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\FullyPreprocessed\Results'
fid = [fpath fname];
load('SegmentedNucleiSegmentationInfo.mat', 'SegmentationResults')
Nucleidata = SegmentationResults;

% Segdata / SegmentationResults variable loaded in from GUI export

%% Preprocessing
% Pulling out the points from the GUI output
[P, fileInfo] = GetPointDetections(Nucleidata);

% Getting the NN distance
Rstats   = GetNuclearDist(P);
D        = {Rstats.Dist2NN};

% Average of nearest
DD       = cellfun(@(x) mean(x), D);
DD2      = reshape(DD, [80 9]);

% Derivative for reshaped array
dDD2     = convn(DD2, [1 0 -1]/2, 'same');

% Where to look for approximate r from NN
avg_dDD2 = mean(dDD2, 2, 'omitmissing');
sqD      = avg_dDD2.^2;

% The minima of sqD is where we can find consistent sizing - dD/dZ ~= 0
[mn, idx] = min(sqD);
Dmin      = DD2(idx, :);
r         = mean(Dmin)/2; % approximate radius in pixels for DBSCAN

% Radius in microns:
% r*.26

%% Clustering Points


% ASK USER FOR METADATA



% Confidence in detection suggests a percentage of the nucleus that needs
% to have been detected to count as a valid detection
CID = 0.5;

% Extracting condition group names from file names
fnames   = {}; 
numfiles = max(P(:,4));
for i = 1:numfiles
    idx = ((i-1)*80)+1; 
    [~, fn] = fileparts(fileInfo{idx});  
    % fnames{i} = fileInfo{idx};
    cyr       = year(datetime);
    idx2      = regexp(fn, ['_' num2str(cyr) '_']);
    idx2      = idx2-1;
    fnames{i} = fn(1:idx2); 
end

% Information on detections after performing 3D clustering with DBSCAN
Nuclei = struct('FilePath', [],...
                'OriginalPoints', [], ...
                'Centroids', [],...
                'CorePoints', [], ...
                'NegatedPoints', [],...
                'PercentExcluded', [],...
                'Cluster', []);

% rescale factor for r - dx/dz gives us a 
rsfactor = .5/.26;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE:
% Need to ensure scaling along z. Incorrect detections (over detection)
% previously occured because points were treated as further away than they
% really are, resulting in insufficient clustering/grouping together. We
% should see increased compaction / a reduction in number of clusters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Looping over files
for f = 1:numfiles
    idx = P(:,4) == f;
    row = find(idx);
    pts = P(row, 1:3);

    % Isotropically scales z to enable proper clustering
    pts(:, 3) = pts(:, 3)/rsfactor;

    % Timing information
    t = tic;

    % Clusters the nuclei in 3D
    [P3D,...
     corepts, ...
     negpts, ...
     discardpercentage,...
     clustergroup] = ClusterNucleiCentroids(pts, r);
    tt = toc(t);

    msg = ['Time to cluster: ' num2str(tt, '%.2f') ' s'];
    disp(msg)

    % Saving information in data struct
    Nuclei(f).FilePath        = fnames{f};
    Nuclei(f).OriginalPoints  = pts;
    Nuclei(f).Centroids       = P3D;
    Nuclei(f).CorePoints      = corepts;
    Nuclei(f).NegatedPoints   = negpts;
    Nuclei(f).PercentExcluded = discardpercentage;
    Nuclei(f).Cluster         = clustergroup;

    % Reporting in CL
    msg = ['Done with: ' num2str(f) '/' num2str(numfiles)];
    disp(msg)
end


%% Saving results of nuclear detection
save NuclearDetections.mat Nuclei r dx dy dz -v7.3 -nocompression


%% Appending results for visualization


% Below is NOT RELEVANT for actual analysis

% AllPoints    = [];
% AllPointsCat = [];
% for i = 1:numfiles
%     pts          = Nuclei(i).Centroids;
%     AllPoints    = vertcat(AllPoints, pts);
%     y            = pts(:,2);
%     x            = pts(:,1);
%     z            = pts(:,3);
%     catshift     = (2048*(i-1));
%     x            = catshift + x;
%     pts2         = [x y z];
%     AllPointsCat = vertcat(AllPointsCat, pts2);
% end
% 
% 
% %% Plotting results from 3D clustering
% 
% 
% % Extracting file names
% fnames   = {}; 
% numfiles = max(P(:,4));
% for i = 1:numfiles
%     idx       = ((i-1)*80)+1; 
%     [~, fn]   = fileparts(fileInfo{idx});  
%     fnames{i} = fileInfo{idx};
% end
% 
% 
% figure

% Max projections and mean projections
% MIPs    = {};
% MeanIPS = {};
% for i = 1:9
%     I            = tiffreadVolume(fnames{i});
%     AllImages{i} = I;
%     MIPs{i}      = max(I, [], 3);
%     MeanIPS      = mean(I, 3);
% end

%% Visualizing all points on images
% if MIPs are computed
% figure
% if ~isempty(AllMIPs)
%     AllMIPs = horzcat(MIPs{:})';
%     imagesc(AllMIPs)
%     colormap gray
%     axis equal
%     hold on
%     plot3(AllPointsCat(:,2), AllPointsCat(:,1), AllPointsCat(:,3), '.r')
% 
%     % Adjusting contrast and lims
%     ax      = gca; 
%     ax.CLim = [0 1500];
%     ax.ZLim = [0 80];
% end


%% Validation of Registration via Phase Correlation

% I = tiffreadVolume('D:/MichaelPalmer_Gillis/2024_04_26 - ATDC5 Cells/IlluminationProfileCorrected/Control1_2024_Apr_26_11.02.45.279_405_0.20502ms_1_Well1BeamCorrected.tif');
% J = tiffreadVolume('D:/MichaelPalmer_Gillis/2024_04_26 - ATDC5 Cells/IlluminationProfileCorrected/Control1_2024_Apr_26_11.08.47.026_405_0.20502ms_1_Well1_BeamCorrected.tif');
% J2 = circshift(J, 8);

% Comparing similarity between slices in the stacks
% numslices         = size(I, 3);
% priorsimilarity   = zeros(numslices, 1);
% postregsimilarity = zeros(numslices, 1);
% for i = 1:numslices
%     im1                  = I(:,:,i);
%     im2                  = J(:,:,i);
%     im3                  = J2(:,:,i);
%     priorsimilarity(i)   = ssim(im1, im2);
%     postregsimilarity(i) = ssim(im1, im3);
%     pct                  = 100*i/numslices;
%     msg                  = ['Percent done: ' num2str(pct, '%.2f')];
%     disp(msg)
% end

%% Plotting Validation of Registration

% figure; plot(priorsimilarity); hold on; plot(postregsimilarity);
% ax = gca;
% ax.YLim = [0.95 1.05];
% yline(1)
% title('Similarity Between Stacks Across Z Improves Post Registration')
% ylabel('Image Simiarlity (SSIM)'); xlabel('Image Slice (z index)');
% legend({'Before Registration', 'Phase Correlation Registration', 'Perfect Similarity (SSIM = 1)'})
% ax.FontSize = 12;




%% Add the path with results

addpath('D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\IlluminationProfileCorrected\Results')

%% Select nuclear channel detections - POST DBSCAN 
load NuclearDetections.mat

%% Select channel of interest segmentation
load Segmentation_488Info.mat SegmentationResults
Ch1 = cat(2, SegmentationResults.SegmentationInfo);

load Segmentation_640Info.mat SegmentationResults
Ch2 = cat(2, SegmentationResults.SegmentationInfo);

% Clean up memory
clear SegmentationResults


%% Getting all file paths

FIDs = [{Ch1(1,:).FilePath} {Ch2(1,:).FilePath}];

% Finding file with z shift (missing z slices)
imds      = imageDatastore(FIDs);
FilesInfo = cellfun(@dir, imds.Files, 'UniformOutput', false);
Filesz    = cellfun(@(x) x.bytes, FilesInfo);

% Largest file - greatest indices
[~, fmx]  = max(Filesz);
fileInfo  = imfinfo(imds.Files{fmx});
Zmx       = numel(fileInfo);
Mmx       = fileInfo(1).Height;
Nmx       = fileInfo(1).Width;

% Smallest file - lesser indices
[~, fmn]  = min(Filesz);
fileInfo  = imfinfo(imds.Files{fmn});
Zmn       = numel(fileInfo);
Mmn       = fileInfo(1).Height;
Nmn       = fileInfo(1).Width;

% Indices to use for region cropping
y = [50 Mmn-50];  % The line scan has some artifacts along z at the ends
x = [1 Nmn];
z = [1 Zmn];

% y = [1 Mmn];




%% Classification of Ch1

% Perform analysis on images
cropregion  = [x' y' z'];
Results_Ch1 = ClassifyExpression(Nuclei, Ch1, r, dx, dz, cropregion, []);

%% Classification of Ch2

% Perform analysis on images
Results_Ch2    = ClassifyExpression(Nuclei, Ch2, r, dx, dz, cropregion, []);

%% Condition names

CNames = {'Control 1', 'Control 2', 'Control 3', ...
          'Day 0 - Well 1', 'Day 0 - Well 2','Day 0 - Well 3',...
          'Day 2 - Well 1', 'Day 2 - Well 2','Day 2 - Well 3'};

for i = 1:9
    Ch1R(i).Expression = Results_Ch1{i}.Expression;
    Ch2R(i).Expression = Results_Ch2{i}.Expression;
    Ch1R(i).F          = Results_Ch1{i}.SumF;
    Ch2R(i).F          = Results_Ch2{i}.SumF;
    Ch1R(i).Fcell      = Results_Ch1{i}.FPerCellA;
    Ch2R(i).Fcell      = Results_Ch2{i}.FPerCellA;
    Ch1R(i).Fvol       = Results_Ch1{i}.FPerCellV;
    Ch2R(i).Fvol       = Results_Ch2{i}.FPerCellV;
    Ch1R(i).Fvox       = Results_Ch1{i}.FPerVox;
    Ch2R(i).Fvox       = Results_Ch2{i}.FPerVox;
end

%%

for i = 1:9
    AVGEXP(i, 1) = mean([Ch1R(i).Expression{:}]);
    AVGVAR(i, 1) = var([Ch1R(i).Expression{:}]);
    AVGEXP(i, 2) = mean([Ch2R(i).Expression{:}]);
    AVGVAR(i, 2) = var([Ch2R(i).Expression{:}]);
end

%% All
ALLRESULTS = [];

for i = 1:9
    Ch1Class           = logical([Ch1R(i).Expression{:}]');
    Ch2Class           = logical([Ch2R(i).Expression{:}]');
    ALLRESULTS(i).CH1  = Ch1Class & ~Ch2Class;
    ALLRESULTS(i).CH2  = Ch2Class & ~Ch1Class;
    ALLRESULTS(i).Both = Ch1Class & Ch2Class;
    ALLRESULTS(i).None = ~Ch1Class & ~Ch2Class;
end


% Creating plots

for i = 1:9
    ALLRESULTS2(i).Ch1 = sum(ALLRESULTS(i).CH1); 
    ALLRESULTS2(i).Ch2 = sum(ALLRESULTS(i).CH2); 
    ALLRESULTS2(i).Both = sum(ALLRESULTS(i).Both); 
    ALLRESULTS2(i).None = sum(ALLRESULTS(i).None);
end

CTypes = {'Col2a Only', 'ColX Only', 'Both', 'Neither'};


for i = 1:9
    NucCount(i) = numel(ALLRESULTS(i).CH1);
end


figure
CH1COUNTS  = [ALLRESULTS2(:).Ch1];
CH2COUNTS  = [ALLRESULTS2(:).Ch2];
BOTHCOUNTS = [ALLRESULTS2(:).Both];
NoneCOUNTS = [ALLRESULTS2(:).None];

Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';
Totals = Totals./NucCount';
Totals = Totals*100;
bar(Totals, 'stacked')
title('Protein Expression')
ylabel('Fraction of Nuclei (%)')

xlabel('Group')
ax = gca;
ax.FontSize = 14;

legend(CTypes)
xticklabels(CNames)



figure
Totals = [CH1COUNTS; CH2COUNTS; BOTHCOUNTS; NoneCOUNTS]';

bar(Totals, 'stacked')
title('Nuclei')
ylabel('Number of Nuclei (Counts)')

xlabel('Group')
ax = gca;
ax.FontSize = 14;

legend(CTypes)
xticklabels(CNames)