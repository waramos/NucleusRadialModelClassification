%% Add the path with results
% 
% addpath('D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\IlluminationProfileCorrected\Results')
% addpath(genpath(cd))
% 
% % Select nuclear channel detections - POST DBSCAN 
% load NuclearDetections.mat
% 
% % Load paths to the files to be loaded in
% load ChannelFilePaths.mat FIDs
% 
% %% 
% % Finding file with z shift (missing z slices)
% imds      = imageDatastore(FIDs);
% FilesInfo = cellfun(@dir, imds.Files, 'UniformOutput', false);
% Filesz    = cellfun(@(x) x.bytes, FilesInfo);
% 
% % Largest file - greatest indices
% [~, fmx]  = max(Filesz);
% fileInfo  = imfinfo(imds.Files{fmx});
% Zmx       = numel(fileInfo);
% Mmx       = fileInfo(1).Height;
% Nmx       = fileInfo(1).Width;
% 
% % Smallest file - lesser indices
% [~, fmn]  = min(Filesz);
% fileInfo  = imfinfo(imds.Files{fmn});
% Zmn       = numel(fileInfo);
% Mmn       = fileInfo(1).Height;
% Nmn       = fileInfo(1).Width;
% 
% %% Cropping
% % Indices to use in lazy loading
% yr = 50;             % Number of y lines to remove
% y  = [yr+1 Mmn-yr];  % The line scan has some artifacts along z at the ends
% x  = [1 Nmn];
% z  = [1 Zmn];
% 
% % Total image sizes
% M = y(2) - y(1) + 1;
% N = x(2) - x(1) + 1;
% Z = z(2) - z(1) + 1;
% 
% % Adjust values for crop region
% y          = [1, M];
% x          = [1, N];
% z          = [1, Z];
% cropregion = [x' y' z'];
% 
% 
% %% Loading in all data
% 
% nfiles = 9;
% nchs   = 2;
% I      = zeros(M, N, Z, nfiles, nchs, 'uint16');
% idx    = 1;
% 
% ntot   = nfiles * nchs;
% 
% for ch = 1:nchs
%     for f = 1:nfiles
%         fid               = imds.Files{idx};
%         J                 = tiffreadVolume(fid, 'PixelRegion', {[y(1) y(2)] [x(1) x(2)] [z(1) z(2)]});
%         I(:, :, :, f, ch) = J;
% 
%         % Verbose reporting on file loading
%         msg = ['Read in: ' newline...
%                fid...
%                newline...
%                'File: ' num2str(idx) '/' num2str(ntot)];
%         disp(msg)
%         idx = idx + 1;
%     end
% end

%% Perform segmentation across different values

% Hyper params of grid search
numsteps  = 5;
stepsz    = -10;

% Parameters 1 and 2: Thresholds for segmentation
p1       = 100;
p2       = 85;
ipvalues = [p1 p2];

% Creates grid of values needed to test sensitivity
Params         = GridSearchVariables(ipvalues, stepsz, numsteps);
[M, N]         = size(Params, [1 2]);
numiterations  = M*N;
iterationcount = 0;

% Min number of voxels required to be within overlap is known
voxthresh = 49;

% Chooses save location for sensitivity analysis results via UI dialog box
% [fname_target, fpath_target] = uiputfile('*.mat', 'Select Save Location');
% [~, fname_target, fext_target] = fileparts(fname_target);


% Loops through the permutations of parameters - param. 1/2 index
for p1idx = 1:1
    for p2idx = 1:1
        % Timing
        t = tic;

        % Iteration/trial # counting
        iterationcount = iterationcount + 1;
        pct            = iterationcount/numiterations;
        pct            = pct*100;

        % Reporting on new iteration across grid and progress
        v1  = Params(1, p1idx, 1);
        v2  = Params(p2idx, 1, 2);
        msg = ['Var 1: ' num2str(v1) newline...
               'Var 2: ' num2str(v2) newline];
        msg = ['New Parameter combination Experiment' newline...
               'Trial: ' num2str(iterationcount)...
               ' of ' num2str(numiterations) newline msg];
        disp(msg)

        % Vector of parameters
        params = [v1 v2];

        % Loop to analyze all files and both channels
        for f = 1:nfiles
            NewResults = [];
            for ch = 1:2
                % Pulling image to analyze
                Im     = I(:, :, :, f, ch);
    
                % Creating mask
                thresh = params(ch);                % Threshold param depends on channel of interest
                Mask   = HardThreshold(Im, thresh); % Segmentation
                Mask   = RefineMask(Mask, 3);       % Cleans up segmentation
    
                % Perform analysis on images
                NewResults = ClassifyExpression(Nuclei(:,f), Mask, voxthresh, r, dx, dz, cropregion, Im);
                Results    = cat(2, Results, NewResults);
            end
        end

        % Appending info
        % if memory is not too intense
        ExpressionResults(p1idx, p2idx).Results  = Results;
        ExpressionResults(p1idx, p2idx).Channel  = ch;
        ExpressionResults(p1idx, p2idx).Trial    = iterationcount;
        ExpressionResults(p1idx, p2idx).Param1   = v1;
        ExpressionResults(p1idx, p2idx).Param2   = v2;


        % Constructing the file path for file to be exported
        fsuffix   = ['_P1_' num2str(v1) '_P2_' num2str(v2)...
                     '_results_' char(datetime('today'))];
        targetFID = [fpath_target fname_target fsuffix fext_target];

        % Progress report
        msg = ['All files and channels analyzed...' newline];
        msg = [msg 'Percent progress: ' num2str(pct, '%.2f') '%'];
        disp(msg)

        % Saving update
        msg = ['Saving results file in: ' newline targetFID '...'];
        disp(msg)

        % Saving
        save(targetFID, 'ExpressionResults', '-v7.3', '-nocompression')

        % Saving update
        msg = 'Results saved!';
        msg = [msg newline...
               'Time for trial: ' toc(t) 's'];
        disp(msg)
    end
end