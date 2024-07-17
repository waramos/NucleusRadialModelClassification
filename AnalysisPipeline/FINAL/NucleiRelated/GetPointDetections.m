function [P, fileInfo] = GetPointDetections(SegmentationResults)
% GETPOINTSDETECTIONS will convert segmentation results into an MxD array
% where M is the number of detections and D is the number of dimensions the
% point detections have (e.g. x, y, z coordinates typical for image stack). 

    % If input is a file path, the data is first loaded
    if (isstring(SegmentationResults) || ischar(SegmentationResults)) && isfile(SegmentationResults)
        load(SegmentationResults, 'SegmentationResults')
    end


    % Pulling out the results
    if isstruct(SegmentationResults) && isfield(SegmentationResults, 'Date') && isfield(SegmentationResults, 'Time')
        Segdata = [SegmentationResults.SegmentationInfo];
    elseif isstruct(SegmentationResults) && isfield(SegmentationResults, 'Results')
        Segdata = SegmentationResults;
    end
    
    % Getting number of z slices and timepoints from the results
    numslices = size(Segdata, 1);
    numframes = size(Segdata, 2);
    numimages = numslices*numframes;
    
    % Reorganizing the data into a single array
    P   = cell(numimages, 2);
    idx = 0;
    for  t = 1:numframes
        for z = 1:numslices

            % New index
            idx    = idx + 1;

            % Always grab file name
            P{idx, 2} = Segdata(z, t).FilePath;

            % adding additional info
            pts    = Segdata(z, t).Results;
            if isempty(pts)
                % In case of a lack of points
                continue
            end
            npts   = size(pts, 1);
            zcoord = repmat(z, [npts, 1]);
            tcoord = repmat(t, [npts, 1]);
            pts    = [pts zcoord tcoord];

            % Filling cell with data
            P{idx, 1} = pts;
            
        end
    end

    % Keeps information separated into points and the corresponding file
    % information if user desires two outputs
    if nargout == 2
        fileInfo = P(:,2);
        P        = cat(1,P{:, 1});
    end
end