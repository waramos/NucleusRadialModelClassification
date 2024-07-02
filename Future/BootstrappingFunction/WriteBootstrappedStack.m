function WriteBootstrappedStack(fid, fov, nsamples, intershuffle, infocusonly)
% WRITEBOOTSTRAPPEDSTACK will write a bootstrapped version of an image
% stack. If the user does not provide a filepath, the function will prompt
% the user to select file(s). If multiple files are selected, they are
% shuffled together as the function will assume that they are different
% conditions meant to be shuffled for a blind approach. The user can
% alternatively toggle the "intershuffle" of the data off and randomly
% sampled stacks will 

    % Window prompts user to select file(s)
    if nargin < 1 || isempty(fid)
        [fname, fpath] = uigetfile('*.tif','MultiSelect','on');

        % Exits, sets up multiple files, or grabs a single file
        if isnumeric(fname)
            return

        elseif iscell(fname)
            fid      = cellfun(@(fn) fullfile(fpath, fn), fname, 'UniformOutput',false);
            numfiles = numel(fid);

        elseif ischar(fname)
            fid      = {fullfile(fpath, fname)};
            numfiles = 1;

        end
    end

    % Default patch size / Field of View is 512x512 pixels
    if nargin < 2 || isempty(fov)
        fov = [512 512];
    end

    % Patch sizing
    mrows = fov(1);
    ncols = fov(2);

    % Default number of samples is 100 patches
    if nargin < 3 || isempty(nsamples)
        nsamples = 100;
    end

    % By default, the data is shuffled together when multiple files exist
    if nargin < 4 || isempty(intershuffle)
        intershuffle = true;
    end

    % Init array and tiff writer (written by Matthew Parent)
    if intershuffle
        J  = zeros(mrows, ncols, nsamples*numfiles, 'uint16');
    else
        J  = zeros(mrows, ncols, nsamples, 'uint16');
    end

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

        % If data is shuffled, it will be placed into the same stack
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
            [focalidx, ~, ~, fprofile] = EstimateFocalPlaneIndices(I);
            Z    = size(I, 3);
            frac = (Z-focalidx)/Z;
            frac = frac*100;
            msg  = ['Focal planes found: ' num2str(frac, '%.2f')...
                    '% found to be in focus (' num2str(focalidx)...
                    '/' num2str(Z) ' slices)'];
            disp(msg)
            J(:,:,idz) = BootstrapImageVolume(I, mrows, ncols, nsamples, focalidx, fprofile, 1e5); % current resample limit is 1e5 attempts at high SNR
        else
            % All slices in volume are considered
            J(:,:,idz) = BootstrapImageVolume(I, mrows, ncols, nsamples);
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