function Results = ClassifyExpressionFast(PC, Clusters, NucleiMask, Ch, I, r, voxthresh, visnhoods, savenhoods, vfreq)
% CLASSIFYEXPRESSION will classify cells as having a protein expressed or
% not having it expressed. This is dependent on some distance value /
% radius, r, that then suggests a cylindrical model out of consideration of
% anistropy.
%
%
% INPUTS:
%
% PC (Mx3 array) - point cloud from 2D detections
%
% NucleiMask (Logical 3D array)
%
% Ch (Logical 3D array) - segmented channel of interest
%
% I (Image stack) - image stack for respective channel of interest
%
% 
%
%
% OUTPUTS:
%
%

    if nargin < 10 || isempty(vfreq)
        vfreq = 100;
    end

    % Whether or not to save the nuclei's local neighborhoods (ROIs)
    if nargin < 9 || isempty(savenhoods)
        savenhoods = false;
    end

    % Whether or not to visualize the nhood grabs in a figure window
    if nargin < 6 || isempty(visnhoods)
        visnhoods = false;
    end

    % Min. num. of voxels needed in the overlap to count as + expression
    if nargin < 3
        voxthresh = [];
    end

    
    if isempty(voxthresh)
        % Histogram based determination
        voxthresh = DetermineMinVoxThresh(Ch);

        % Verbose reporting
        msg = ['Threshold for positive classification established as: ' ...
               newline...
               num2str(voxthresh) ' Voxels'];
        disp(msg)
    end

    % Init results array
    Results = [];

    % Verbose reporting
    msg = 'Point cloud adjusted';
    disp(msg)

    % Determining number of clusters and their labels
    Labels   = unique(Clusters); 
    nLabels  = numel(Labels);

    % Results per file
    Cluster.PixIdx     = cell(nLabels, 1);
    Cluster.NHood      = cell(nLabels, 1);
    Cluster.Masks      = cell(nLabels, 1);
    Cluster.Expression = cell(nLabels, 1);
    Cluster.SumF       = cell(nLabels, 1);
    Cluster.FPerCellA  = cell(nLabels, 1);
    Cluster.FPerCellV  = cell(nLabels, 1);
    Cluster.FPerVox    = cell(nLabels, 1);

    % Area of cell
    a = pi*(r^2);

    % Counts number of clusters that are thrown out
    badcluster = 0;

    % Sizing information for bound checking the PC to Mask conversion 
    sz = size(Ch, [1 2 3]);

    % Setup visualization figure
    if visnhoods
        noplot       = true;
        fig          = figure;
        fig.Units    = 'Normalized';
        fig.Position = [0.25 0.25 0.5 0.5];
        ax           = axes(fig);
        ax.Units     = 'Normalized';
        ax.Position  = [0.1 0.1 0.8 0.8];
    end


    % Looping through each cluster, i.e. PC of a single cell
    for i = 1:nLabels
        % Index in original PC for given label
        idx       = Clusters == i;

        if nnz(idx) == 0
            badcluster = badcluster + 1;
            msg = ['Bad cluster index: ' num2str(i) newline ...
                   'Number of bad clusters: ' num2str(badcluster) newline];
            disp(msg)
            continue
        end

        % Grab points for local nhood
        [rows, ~] = find(idx);
        pts       = PC(rows, :);

        % Volume is area * number of planes
        numslices = range(pts(:,3));
        vol       = a*numslices;

        % Pulling out local window for each nucleus
        [y, x, z]         = GetLocalCellNHood(pts, r, sz);
        Cluster.PixIdx{i} = {y, x, z};
        Im                = I(y, x, z);
        NucMask           = NucleiMask(y, x, z);
        ChMask            = Ch(y, x, z);

        % Holding image and masks
        if savenhoods
            Cluster.NHood{i} = Im;
            Cluster.Masks{i} = {NucMask ChMask};
        end

        % Overlapping region
        Overlap               = NucMask & ChMask;
        voxcount              = nnz(Overlap);
        Cluster.Expression{i} = voxcount > voxthresh;

        % Fluorescence
        F                    = Im(Overlap);
        F                    = sum(F);
        Cluster.SumF{i}      = F;
        Cluster.FPerCellA{i} = F/a;
        Cluster.FPerCellV{i} = F/vol;
        Cluster.FPerVox{i}   = F/voxcount;

        
        if (mod(i, vfreq) == 0) || i==nLabels
            % Verbose reporting
            pct = 100*(i/nLabels);
            msg = ['Cells processed: ' num2str(i) '/' num2str(nLabels)...
                   newline num2str(pct, '%.2f') '% Done' newline];
            disp(msg)

            % Visualization of nhood's maximum projection
            if visnhoods
                Maxp = max(Im, [], 3);
    
                % Checks if plot was created before
                if noplot
                    iH = imagesc(ax, Maxp); 
                    axis image
                    colormap bone
                else
                    iH.CData = Maxp;
                end
    
                % Title update to the axes
                t = ['Percent Done: ' num2str(pct, '%.2f')];
                title(t)
                drawnow
            end
        end  
    end

    % Concatenate results vertically in the data struct
    Results = cat(1, Results, Cluster);
end