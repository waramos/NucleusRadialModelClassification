%% Load in a subset of the files


% Starting and ending indices along image dimensions
y1 = 1;
y2 = 2048;
x1 = 1;
x2 = 2048;
z1 = 1;
z2 = inf;

% Neigborhood to load in for registration
rows = [y1 y2];
cols = [x1 x2];
z    = [1 inf];

% Selecting the folder of images to process and register
selpath = uigetdir(cd, 'Select folder to process');
imds    = imageDatastore(selpath, "FileExtensions",'.tif', 'ReadFcn', @(x) ReadSubSetOfStack(x, rows, cols));


%% User selection of registration channel

% Extracting file names to determine which channels were imaged
fnames         = imds.Files;
strexp         = '_[0-9]{3,3}_'; % three numbers used to describe the channel
[s_idx, e_idx] = regexp(fnames, strexp);

% Number of files
nfiles = numel(fnames);
ch     = {};
for i = 1:nfiles
    % Extracting information regarding channels
    fid    = fnames{i};
    chused = fid(s_idx{i}:e_idx{i});
    ch{i}  = chused(2:end-1);        % eliminates underscores
end

% All unique channels imaged
ch = unique(ch);

% User to select the reference channel
pmsg  = 'Select reference channel for alignment of other data.';
chidx = listdlg('ListString', ch, ...
                'PromptString', pmsg,...
                'SelectionMode', 'single', ...
                'InitialValue', 1);
refch = ch{chidx};

%% Registration of Images

% Splitting data into the reference files and files to be registered
isrefch   = contains(fnames, refch);
reffiles  = fnames(isrefch);
files2reg = fnames(~isrefch);

% Initialze regitration class
REG = Registration;

% Initialize the tiff writing class
FW = FWTiff;

% Loop will read in groups of two for registration
ngroups = numel(reffiles)/2;
shifts  = zeros(ngroups, 3);
rdims   = [1 2 3];
for i = 1:ngroups
    % File index
    fidx  = 2*(i-1) + 1;
    fidx2 = fidx + 1;

    % Extract the file IDs
    fid1 = reffiles{fidx};
    fid2 = reffiles(fidx2);

    % Read in the image volumes.
    I1 = ReadSubSetOfStack(fid1, rows, cols);
    I2 = ReadSubSetOfStack(fid2, rows, cols);

    % Progress reporting
    msg = ['Data Loaded'...
           newline...
           'File group: ' num2str(i) '/' num2str(ngroups)];
    disp(msg)

    % % Finding shifts from cross power spectrum
    [r, Q]       = REG.FourierShiftND(I1, I2, rdims);
    shifts(i, :) = -r;
    r            = shifts(i, :);

    % Progress reporting
    msg = ['Registration calculated'...
           newline...
           'Shifts found : ' num2str(shifts(i, :)) ' along dimensions ' num2str(rdims)];
    disp(msg)


    % Only need to shift the second file relative to the first
    fid = files2reg{fidx2};
    I   = tiffreadVolume(fid);
    
    % Checks for certain shift that might not normally be expected
    haszshift = r(3)~=0;
    hasxshift = r(2)~=0;

    % Will force x shift to be 0 ?
    if hasxshift
        xshift_og = shifts(i, 2);
        r(2)      = 0;
        msg       = ['Original x shift: ' num2str(xshift_og)];
        disp(msg)
    end

    % Shifts the image accordingly
    I   = circshift(I, r);

    % Will discard some slices
    if haszshift
        if r(3) > 0
            z1 = r(3)+1;
            z2 = size(I, 3);
        else
            z1 = 1;
            z2 = size(I, 3) + r(3);
        end

        % Adjusts the z planes of the image
        I   = I(:,:,z1:z2);
        msg = ['Z shift was found and some slices were removed.' ...
               newline...
               'Z Shift: ' num2str(r(3))...
               newline...
               'z slices: ' num2str(z1) ' to ' num2str(z2)];
        disp(msg)
    end

    % Creating new file ID
    [fpath, fname, fext] = fileparts(fid);
    newpath = [fpath filesep 'Registered'];
    newfid  = [newpath filesep fname '_Registered' fext];

    % In case the new path does not yet exist
    if ~exist(newpath, "dir") == 7
        mkdir(newpath)
    end

    % Writing new file under the desired FID
    FW.Open(newfid)
    FW.Write(I)
    FW.Close
    
    % Minor update
    pct = i/ngroups;
    pct = 100*pct;
    pct = num2str(pct);
    msg = ['File written' ...
           'Progress: ' pct];
    disp(msg)
    
end




%% Functions
function I = ReadSubSetOfStack(fid, y, x, z)
    if nargin < 4
        z = [1 inf];
    end
    if nargin < 3
        x = [1 inf];
    end
    if nargin < 2
        y = [1 inf];
    end

    % Reads in the subset of the image stack
    I = tiffreadVolume(fid, "PixelRegion", {y, x, z});
end