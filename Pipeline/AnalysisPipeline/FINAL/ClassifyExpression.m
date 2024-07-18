function Results = ClassifyExpression(NInfo, ChInfo, r, dx, dz, cropregion, I, visnhoods, savenhoods)
% CLASSIFYEXPRESSION will classify cells as having a protein expressed or
% not having it expressed. This is dependent on some distance value /
% radius, r, that then suggests a cylindrical model out of consideration of
% anistropy.


    if nargin < 9 || isempty(savenhoods)
        savenhoods = false;
    end

    if nargin < 8 || isempty(visnhoods)
        visnhoods = false;
    end

    if nargin < 7 || isempty(I) 
        readfilesin = true;
    else
        readfilesin = false;
    end

    % isotropic scaling
    if dx ~= dz
        rsfactor = dz/dx;
    else
        rsfactor = 1;
    end

    % Init results
    Results = [];

    % Check if channel info is a struct or an array (mask stack)
    isstack = ~isstruct(ChInfo) && isnumeric(ChInfo);

    % Number of files to analyze
    nFiles = numel(NInfo);
    for f = 1:nFiles
        % Verbose reporting
        msg = ['Analyzing file: ' num2str(f) '/' num2str(nFiles)];
        disp(msg)

        if readfilesin && ~isstack
            % Reading in the original image file
            I   = tiffreadVolume(ChInfo(1,f).FilePath);
            msg = 'Loading Channel data in';
            disp(msg)
        end

        % Size information
        sz = size(I, [1 2 3]);

        % Default "crop" grabs everything
        if nargin < 6 || isempty(cropregion)
            cropregion = [ones(1,3); sz];
        end

        % Pulling out original PC and the cluster labels
        PC       = NInfo(1, f).OriginalPoints;
        Clusters = NInfo(1, f).Cluster;

        % Corrections and cropping as needed - puts PC in image coords
        PC             = AdjustPointCloud(PC, sz, rsfactor);
        [PC, Clusters] = CropPC(PC, Clusters, cropregion);

        % Binary image from nuclei detections
        Nuclei = Points2Cylinders(PC, r, sz);
        Nuclei = CropVol(Nuclei, cropregion);

        % Binary image from channel segmentation
        if ~isstack
            Ch = cat(3, ChInfo(:, f).Results);
        else
            Ch = ChInfo;
            clear ChInfo
        end

        % Cropping as needed
        Ch = CropVol(Ch, cropregion);

        % Verbose reporting
        msg = 'Determining protein expression threshold';
        disp(msg)

        % Minimum number of voxels needed in the overlap to count as
        % positive classification
        thresh = DetermineMinVoxThresh(Ch);

        % Verbose reporting
        msg = ['Threshold for positive classification established as: ' ...
               newline...
               num2str(thresh) ' Voxels'];
        disp(msg)

        

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

        % Visualizes neighorhoods if desired
        if visnhoods
            ff     = figure;
            ax     = axes(ff);
            noplot = true;
        end        

        % Sizing information for bound checking the PC to Mask conversion 
        sz2 = size(Ch, [1 2 3]);

        for i = 1:nLabels
            % Index in original PC for given label
            idx       = Clusters == i;

            if nnz(idx) == 0
                badcluster = badcluster + 1;
                msg = ['Bad cluster index: ' num2str(i) newline ...
                       'Number of bad clusters: ' num2str(badcluster)];
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
            [y, x, z]         = GetLocalCellNHood(pts, r, sz2);
            Cluster.PixIdx{i} = {y, x, z};
            Im                = I(y, x, z);
            NucMask           = Nuclei(y, x, z);
            ChMask            = Ch(y, x, z);

            % Holding image and masks
            if savenhoods
                Cluster.NHood{i} = Im;
                Cluster.Masks{i} = {NucMask ChMask};
            end

            % Overlapping region
            Overlap   = NucMask & ChMask;
            voxcount  = nnz(Overlap);
            Cluster.Expression{i} = voxcount > thresh;

            % Fluorescence
            F                    = Im(Overlap);
            F                    = sum(F);
            Cluster.SumF{i}      = F;
            Cluster.FPerCellA{i} = F/a;
            Cluster.FPerCellV{i} = F/vol;
            Cluster.FPerVox{i}   = F/voxcount;

            % Verbose reporting
            pct = 100*(i/nLabels);
            msg = ['Percent nuclei done: ' num2str(pct, '%.2f')];
            disp(msg)

            % Visualization of maximum projection every 10 iterations
            if mod(i, 10) == 0 && visnhoods
                Maxp = max(Im, [], 3);

                % Checks if plot was created before
                if noplot
                    iH = imagesc(ax, Maxp); 
                    axis image; 
                    colormap bone
                    noplot = false;
                else
                    iH.CData = Maxp;
                end

                % Title update to the axes
                t = ['Percent Done: ' num2str(pct, '%.2f')];
                title(t)
                drawnow
            end
            
        end

        % Consolidate results
        Results = cat(1, Results, Cluster);

        msg = ['Done with file: ' num2str(f) '/' num2str(nFiles)];
        disp(msg)

    end

end



function PC = AdjustPointCloud(PC, sz, rsfactor)
% ADJUSTPOINTCLOUD ensures the point cloud is rounded to the nearest pixel,
% then corrected by the rescale factor, and lastly bound checks the final
% point cloud to make sure it can index the image accurately.
    % Converting size vector to x y z ordering
    sz      = [sz(2) sz(1) sz(3)];
    PC(:,3) = PC(:,3) * rsfactor;
    PC      = round(PC);
    PC      = max(1, PC);
    PC      = min(PC, sz);
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

    % Adjust the point clouds by the shift
    PC = PC - [x1 y1 z1];
end



function [y, x, z, sz] = GetLocalCellNHood(pts, r, sz)
% GETLOCALCELLNHOOD constructs a local neighborhood based on a pointcloud
% of detections
    [x1, x2] = bounds(pts(:,1));
    [y1, y2] = bounds(pts(:,2));
    [z1, z2] = bounds(pts(:,3));

    % X Y Z coordinates for the local neighborhood/window to grab
    y = (y1-r):(y2+r);
    x = (x1-r):(x2+r);
    z = (z1-r):(z2+r);
    y = round(y);
    x = round(x);
    z = round(z);

    % Bound check
    y = min(y, sz(1));
    x = min(x, sz(2));
    z = min(z, sz(3));
    y = max(y, 1);
    x = max(x, 1);
    z = max(z, 1);

    try
        dy = y(end)-y(1);
        dx = x(end)-x(1);
        dz = z(end)-z(1);
    catch
        keyboard
        disp('')
    end

    % Sizing information
    sz = [dy dx dz];
end


function Nuclei = Points2Cylinders(pts, r, sz)

    % Round to nearest pixel - new var enables comparison
    PC_pix = round(pts);
    x      = PC_pix(:,1);
    y      = PC_pix(:,2);
    z      = PC_pix(:,3);

    % Convert point cloud to linear indices corresponding to pixels
    idx = sub2ind(sz, y, x, z);

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


function thresh = DetermineMinVoxThresh(Mask)
% DETERMINEMINVOXTHRESH determines the minimum number of voxels needed to
% be overlapping for a positive classification.
    CC         = bwconncomp(Mask, 8);
    RP         = regionprops(CC, 'Area');
    A          = [RP.Area];
    [N, edges] = histcounts(A, 1:1000);
    [~, idx]   = max(N);
    thresh     = edges(idx);
end

function isabove = AboveThresh(O, thresh)
    isabove = nnz(O)>thresh;
end
