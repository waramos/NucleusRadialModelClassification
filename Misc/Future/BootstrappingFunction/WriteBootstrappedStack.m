function Indices = WriteBootstrappedStack(fid, fov, nsamples, intershuffle, infocusonly)
% WRITEBOOTSTRAPPEDSTACK will write a bootstrapped version of an image
% stack. If the user does not provide a filepath, the function will prompt
% the user to select file(s). If multiple files are selected, they are
% shuffled together as the function will assume that they are different
% conditions meant to be shuffled for a blind approach. The user can
% alternatively toggle the "intershuffle" of the data off and randomly
% sampled stacks will instead by 

    % Window prompts user to select file(s)
    if nargin < 1 || isempty(fid)
        [fname, fpath] = uigetfile('*.tif','MultiSelect','on');

        % Exits, sets up multiple files, or grabs a single file
        if isnumeric(fname)
            disp('No files selected')
            return
        else
            fid = fname;
        end
    end

    % In case multiple files are selected or passed in as a cell array
    if iscell(fid)
        fid      = cellfun(@(fn) fullfile(fpath, fn), fid, 'UniformOutput',false);
        numfiles = numel(fid);

    elseif ischar(fid)
        fid      = {fullfile(fpath, fid)};
        numfiles = 1;
    end

    if nargin < 2 || isempty(fov)
        fov = [64 64 3];
    end

    % Default number of samples is 100 patches
    if nargin < 3 || isempty(nsamples)
        nsamples = 100;
    end

    % By default, the data is shuffled together when multiple files exist
    if nargin < 4 || isempty(intershuffle)
        intershuffle = true;
    end

    % Extract dims
    mrows = fov(1);
    ncols = fov(2);

    % Intershuffle will shuffle files around so user does not know about
    % conditions
    if intershuffle
        J  = zeros(mrows, ncols, nsamples*numfiles, 'uint16');
    else
        J  = zeros(mrows, ncols, nsamples, 'uint16');
    end

    % Class written by Matthew Parent for fast tiff writing
    FW = FWTiff;

    % Loops through all files
    for f = 1:numfiles
        % Load in original data
        FID        = fid{f};
        [fpath,...
         fname,...
         fext]     = fileparts(FID);
        I          = tiffreadVolume(FID);

        % Update on progress in loading
        msg = ['File ' num2str(f) '/' num2str(numfiles) ' loaded.'];
        disp(msg)

        % Checks that init array is same data type
        if f ==1 && ~isa(I, class(J))
            J = cast(J, class(I));
        end

        % If data is shuffled, it will be placed into the same stack and
        % the stack is not written until the very end. Otherwise, each file
        % results in a new stack being written
        if intershuffle
            z1       = (nsamples*(f-1))+1;
            z2       = (f*nsamples);
            idz      = z1:z2;
        else
            idz      = 1:nsamples;
        end

        % Bootstrap data and assign into init. array
        if infocusonly
            % User requests ONLY the data in focus
            [focalidx, ...
             ~, ...
             ~, ...
             fprofile] = EstimateFocalPlaneIndices(I);
            Z        = size(I, 3);
            focalidx = max(focalidx, [1 1]);
            nz       = max(1, focalidx(2)-focalidx(1));
            frac     = (Z-nz)/Z;
            frac     = frac*100;
            msg      = ['Focal planes found: ' num2str(frac, '%.2f')...
                        '% found to be in focus (' num2str(nz)...
                        '/' num2str(Z) ' slices)'];
            disp(msg)
            [J(:,:,idz), Indices(f)] = BootstrapImageVolume(I, fov, nsamples, true, focalidx, fprofile, 1e3); % current resample limit is 1e3 attempts at high SNR
        else
            % All slices in volume are considered
            [J(:,:,idz), Indices(f)] = BootstrapImageVolume(I, fov, nsamples);
        end

        % Writes separate stacks when unshuffled
        if ~intershuffle
            newfname = [fname fext];
            newfid   = fullfile(fpath, newfname);
            FW.Open(newfid)
            FW.Write(J)
            FW.Close

            % Exit function
            return
        end
    end

    % Shuffled data requires user to specify where to save
    if intershuffle
        % New filepath
        [fname, fpath] = uiputfile('*.tif', 'Save location', fpath);
        newfid         = fullfile(fpath, fname);

        % Shuffle the patches
        numz = size(J, 3);
        idz  = randperm(numz);
        J    = J(:,:, idz);

        % Write the new stack
        FW.Open(newfid)
        FW.Write(J)
        FW.Close
    end
end


function [idz, d_df, df, fprofile] = EstimateFocalPlaneIndices(I)
% ESTIMATEFOCALPLANEINDICES will estimate which z slices from an image
% volume are in focus. This is done by first computing the average signal
% across each plane. Then, the gradient of the average fluorescence with
% respect to z is computed. This is smoothed and the point at which the
% greatest gradient increase occurs is approximated to be where content
% first most noticeably comes into focus.
    
% Fluorescence change
    % Average signal profile across the stack
    fprofile = mean(I, [1 2]);
    fprofile = squeeze(fprofile);
    
    % dF/dz
    df = conv(fprofile, [1 0 -1]/2, 'same');
    df = smooth(df(2:end-1));

    % Approx inflection point where things begin to be in focus
    d_df     = conv(df, [1 0 -1]/2, 'same');
    [~, idz1] = max(d_df(2:end-1));
    [~, idz2] = min(d_df(2:end-1));

    % Start/End of focus - adjusted given clipping at ends of df and ddf
    idz = [idz1 idz2] + 2;
end