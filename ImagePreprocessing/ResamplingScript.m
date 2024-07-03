%% Resampling all data


% Folders with data
% Parent directory, ch1 directory, registered images' directory
pdir      = 'D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells';
ch1dir    = [pdir filesep 'IlluminationProfileCorrected'];
regdir    = [ch1dir filesep 'Registered'];
targetdir = [pdir filesep 'FullyPreprocessed'];

% Makes the target directory if it does not already exist
if exist(targetdir, 'dir') ~= 7
    mkdir(targetdir)
end

% Getting filenames from imds
% Reference, Ch1, and Registered files
reffiles = imageDatastore(pdir, 'FileExtensions', '.h5');
ch1files = imageDatastore(ch1dir, 'FileExtensions', '.tif');
regfiles = imageDatastore(regdir, 'FileExtensions', '.tif');

% Extracting the file names
reffiles = reffiles.Files;
ch1files = ch1files.Files;
regfiles = regfiles.Files;


%% Adjusting file names

% Ensure channels are appropriately labeled and we only pull files for a
% given channel
isch1      = contains(ch1files, '_488_');
ch1files   = ch1files(isch1);

isregch    = contains(regfiles, '_640_');
regfiles   = regfiles(isregch);

% Ref files will contain info regarding the voxel size for either channel
iseither = contains(reffiles, '_640_') | contains(reffiles, '_488_');
reffiles = reffiles(iseither);

% Joins both cell arrays to just load all data
files2load = [ch1files; regfiles];


%% Resampling

% Checking the imaging parameters - voxel size info
nfiles    = numel(reffiles);
pixelsize = zeros(nfiles, 1);
zstep     = zeros(nfiles, 1);
for i = 1:nfiles
    % Getting imaging parameters
    reffid       = reffiles{i};
    mdata        = GetH5ImageParams(reffid);
    pixelsize(i) = mdata.pixelsize;
    zstep(i)     = mdata.zstep;
end

% Ensure same voxel size was used, otherwise, comparisons cannot be made
pixelsize = unique(pixelsize);
samedxy   = numel(pixelsize) == 1;
zstep     = unique(zstep);
samedz    = numel(zstep) == 1;

if ~samedxy && ~samedz
    error(['Images have different voxel sizes. ' ...
           'Consider matching imaging parameters.'])
end

% Memory efficiency automatically considered
mem     = memory;
freemem = mem.MemAvailableAllArrays;

% Checks the first file to get an idea of memory needed to interpolate
fileInfo = dir(files2load{1});
fileSize = fileInfo.bytes;

if fileSize < freemem
    memoryeff = false;
else
    memoryeff = true;
end

% Fast Write Tiff class to write resampled files
FW = FWTiff;

% Load in the image stacks
nfiles = numel(files2load); % Number of files to load in
for i = 1:1
    % Load file
    fid    = files2load{i};
    I      = tiffreadVolume(fid);

    % Resampled stack
    J      = IsoResampleImageStack(I, pixelsize, zstep, false);

    % New file path
    [~, fname, ~] = fileparts(fid);
    newfid        = [targetdir fname '.tif'];

    % Writing new file
    FW.Open(newfid)
    FW.Write(J)
    FW.Close
end


