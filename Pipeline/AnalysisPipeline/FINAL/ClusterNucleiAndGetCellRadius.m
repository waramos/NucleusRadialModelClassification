function [Nuclei, r, D] = ClusterNucleiAndGetCellRadius(Nucleidata)
% CLUSTERNUCLEIANDGETCELLRADIUS will pull out a pointcloud representing
% detections of nuclei across 2D image planes from a 3D image stack. It
% will do this in an automated way based off of information regarding known
% imaging parameters from the original files' metadata. This information
% enables scaling the z dimension so that it isotropically matches up with
% x and y. This is then clustered with DBSCAN using the distance, r, found
% from the first nearest neighbor of the point.

    % Pulling out the points from the GUI output
    [P, fileInfo] = GetPointDetections(Nucleidata);
    
    % Getting the NN distance
    Rstats   = GetNuclearDist(P);
    D        = {Rstats.Dist2NN};
    
    % Average (across slices) of dist to nearest neighbor
    DD       = cellfun(@(x) mean(x), D);
    DD2      = reshape(DD, [80 9]);
    
    % dD/dz - change in nearest neighbor distance as z slice changes
    dDD2     = convn(DD2, [1 0 -1]/2, 'same');
    
    % Average the change in NN dist. across files
    avg_dDD2 = mean(dDD2, 2, 'omitmissing');
    sqD      = avg_dDD2.^2;
    
    % Minima of sqD suggests consistent sizing: dD/dZ ~= 0
    [~, idx] = min(sqD);
    Dmin     = DD2(idx, :);
    r        = mean(Dmin)/2; % Approximate cell radius
    
    % Extracting condition group names from file names
    fnames   = {};              % init cell array
    numfiles = max(P(:,4));     % Max file index
    cyr      = year(datetime);  % Current year

    % Looping over the file names
    for i = 1:numfiles
        idx       = ((i-1)*80)+1;
        fid       = fileInfo{idx};
        [~, fn]   = fileparts(fid);  
        FIDs{i}   = fileInfo{idx}; % If full file name desired

        % Only grabs the condition name
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
        Nuclei(f).FilePath        = FIDs{f};
        Nuclei(f).ConditionName   = fnames{f};
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


%% Need to allow user to save to specific place
% Saving results of nuclear detection as well as meta values
% save NuclearDetections.mat Nuclei r dx dy dz -v7.3 -nocompression

end