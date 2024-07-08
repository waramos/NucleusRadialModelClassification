
% Segdata / SegmentationResults variable loaded in from GUI export

% Pulling out the points from the GUI output
% [P, fileInfo] = GetPointDetections(Segdata);

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

%% Clustering Points

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
    idx2 = regexp(fn, '_2024_');
    idx2 = idx2-1;
    fnames{i} = fn(1:idx2); 
end

% Information on detections after performing 3D clustering with DBSCAN
Nuclei2 = struct('FilePath', [],...
                'OriginalPoints', [], ...
                'Centroids', [],...
                'CorePoints', [], ...
                'NegatedPoints', [],...
                'PercentExcluded', [],...
                'Cluster', []);

% dx/dz gives us a rescale factor for r
rsfactor = 0.26/.5;
r_dxdz   = r*rsfactor;

% Looping over files
for f = 1:numfiles
    idx = P(:,4) == f;
    row = find(idx);
    pts = P(row, 1:3);

    t = tic;
    [P3D,...
     corepts, ...
     negpts, ...
     discardpercentage,...
     clustergroup] = ClusterNucleiCentroids(pts, r_dxdz);
    tt = toc(t);

    msg = ['Time to cluster: ' num2str(tt, '%.2f') ' s'];
    disp(msg)

    % Saving information in data struct
    Nuclei2(f).FilePath        = fnames{f};
    Nuclei2(f).OriginalPoints  = pts;
    Nuclei2(f).Centroids       = P3D;
    Nuclei2(f).CorePoints      = corepts;
    Nuclei2(f).NegatedPoints   = negpts;
    Nuclei2(f).PercentExcluded = discardpercentage;
    Nuclei2(f).Cluster         = clustergroup;

    % Reporting in CL
    msg = ['Done with: ' num2str(f) '/' num2str(numfiles)];
    disp(msg)
end

%% Appending results for visualization

AllPoints    = [];
AllPointsCat = [];
for i = 1:numfiles
    pts          = Nuclei2(i).Centroids;
    AllPoints    = vertcat(AllPoints, pts);
    y            = pts(:,2);
    x            = pts(:,1);
    z            = pts(:,3);
    catshift     = (2048*(i-1));
    x            = catshift + x;
    pts2         = [x y z];
    AllPointsCat = vertcat(AllPointsCat, pts2);
end


%% Plotting results from 3D clustering

figure

AllMIPs = horzcat(MIPs{:})';
imagesc(AllMIPs)
colormap gray
axis equal
hold on
plot3(AllPointsCat(:,2), AllPointsCat(:,1), AllPointsCat(:,3), '.r')

% Adjusting contrast and lims
ax      = gca; 
ax.CLim = [0 1500];
ax.ZLim = [0 80];


%% Making spheres
