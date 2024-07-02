function V = LoadAllVolumes(fpaths)
    % Getting all files names
    if isempty(fpaths)
        [fnames, fpath] = uigetfile('*.tif', 'Load', cd, 'MultiSelect', 'on');
        if isnumeric(fpath)
            disp('No files loaded.')
            return
        end
        fpaths          = cellfun(@(fn) fullfile(fpath, fn), fnames, 'UniformOutput', false);
    elseif ~iscell(fpaths)
        fpaths = {fpaths};
    end

    % Metadata revision
    mdata           = imfinfo(fpaths{1});
    
    % Sizing
    M          = mdata(1).Height;
    N          = mdata(1).Width;
    Z          = numel(mdata);
    numfiles   = numel(fpaths);
    V          = zeros(M, N, Z, numfiles, 'uint16');

    % loop for appending the data
    t0 = tic;
    for f = 1:numfiles
        t   = tic;
        fid = fpaths{f};
        V(:,:,:,f) = tiffreadVolume(fid);
    
        % CL Updates
        tt   = toc(t);
        msg  = ['Image ' num2str(f) ' of ' num2str(numfiles) ' loaded and appended.'];
        tmsg = ['Time to load: ' num2str(tt) ' seconds.'];
        msg  = [msg newline tmsg];
        disp(msg)
    end
    
    % Loading done: Final total time update
    t0   = toc(t0);
    tmsg = ['Time to load: ' num2str(t0) ' seconds.'];
    tmsg = [tmsg newline 'Loading data done.'];
    disp(tmsg)
end