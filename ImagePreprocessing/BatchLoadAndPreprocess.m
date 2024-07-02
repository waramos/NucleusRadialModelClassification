% Getting all files names
[fnames, fpath] = uigetfile('*.tif', 'Load', 'D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\', 'MultiSelect', 'on');
if isnumeric(fpath)
    disp('No files loaded.')
    return
end

% Appending the file directory to each file name
fpaths          = cellfun(@(fn) fullfile(fpath, fn), fnames, 'UniformOutput', false);
mdata           = imfinfo(fpaths{1});

% Sizing information
M          = mdata(1).Height;
N          = mdata(1).Width;
Z          = numel(mdata);
numfiles   = numel(fpaths);
AllVolumes = zeros(M, N, Z, numfiles, 'uint16');

%% Batch load

% loop for appending the data
t0 = tic;
for f = 1:numfiles
    t   = tic;
    fid = fpaths{f};
    AllVolumes(:,:,:,f) = tiffreadVolume(fid);

    % CL Updates
    tt   = toc(t);
    msg  = ['Image ' num2str(f) ' of ' num2str(numfiles) ' loaded and appended.'];
    tmsg = ['Time to load: ' num2str(tt) ' seconds.'];
    msg  = [msg newline tmsg];
    disp(msg)
end

% Final total time update
t0   = toc(t0);
tmsg = ['Time to load: ' num2str(t0) ' seconds.'];
tmsg = [tmsg newline 'Loading data done.'];
disp(tmsg)

% Loading done

%% Batch beam correct

% Fast Write Tiff object quickly writes new tif files
FW     = FWTiff;
labels = {'Controls', 'Day0', 'Day2'};
V      = zeros(size(AllVolumes, [1 2 3]), 'uint16');
for i = 1:numfiles
    oldfid = fpaths{i};
    [fpath, fname, fext] = fileparts(oldfid);

    l1 = i<=3;
    l2 = (3<i) && i<=6;
    l3 = i>6;
    trialnum  = mod(i, 3);

    if trialnum == 0
        trialnum = 3;
    end

    fname  = [fname '_Well' num2str(trialnum) '_BeamCorrected'];
    newfid = [fpath filesep fname fext];

    V = AllVolumes(:,:,:,i);
    V = QuickBeamCorrection(V);

    FW.Open(newfid)
    FW.Write(V)
    FW.Close

    prog = 100*i/numfiles;
    prog = ['Percent done: ' num2str(prog)];
    disp(prog)
end