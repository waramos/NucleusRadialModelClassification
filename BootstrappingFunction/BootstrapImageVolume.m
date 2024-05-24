function J = BootstrapImageVolume(I, mrows, ncols, numpatches, idz, fprofile, resamplelimit)
% BOOSTRAPIMAGEVOLUME bootstraps a 3D image stack with some field of view
% based on number of specified rows and columns. It will produce samples as
% patches with the requested field of view and as many samples as are
% requested.

    % Fluorescence profile across the stack for consideration of sampling
    if nargin < 5 || isempty(idz)
        idz = 1;
    end

    if nargin < 6 || isempty(fprofile)
        fprofile = [];
    end

    % Only grab focal planes
    I = I(:,:,idz:end);

    % Minimum window sizing
    mn = 64;

    % Get image bounds
    [M, N, Z] = size(I, [1 2 3]);

    % In case percentage of dimension requested
    if mrows < 1
        mrows = mrows*M;
    end
    if ncols < 1
        ncols = ncols*N;
    end

    % Ensure integer valued FOV sizing
    mrows = round(mrows);
    ncols = round(ncols);

    % Bounds
    mrows = max(mrows, mn);
    ncols = max(ncols, mn);
    newM  = M - mrows;
    newN  = N - ncols;

    % Sampling indices
    idy = randi(newM, [1 numpatches]);
    idx = randi(newN, [1 numpatches]);
    idz = randi(Z, [1 numpatches]);

    % Maximum number of resamples allowed when user only wants patches
    % with content in focus
    if nargin < 6 || isempty(resamplelimit) 
        resamplelimit = 1000;
    end

    % Init shuffled image array
    ogclass = class(I);
    J       = zeros(mrows, ncols, numpatches, ogclass);
    for i = 1:numpatches
        % Init patch resampling params each slice
        goodpatch     = false;
        resamplecount = 0;

        while ~goodpatch && resamplecount < resamplelimit
            % Resample count update
            resamplecount = resamplecount + 1;
            if resamplecount > 1 && mod(resamplecount, 100) == 0
                msg = ['Low SNR patch. Resampling...' ...
                       newline ...
                       'Resample Iteration: ' num2str(resamplecount) ...
                       newline...
                       'Limit: ' num2str(resamplelimit)];
                disp(msg)
            end

            % Start and end indices
            y1 = idy(i);
            x1 = idx(i);
            y2 = y1 + mrows - 1;
            x2 = x1 + ncols - 1;
    
            % Patch indices
            y  = y1:y2;
            x  = x1:x2;
            z  = idz(i);
    
            % Grabbing image patch;
            impatch = I(y, x, z);
            if ~isempty(fprofile)
                goodpatch = CheckSNR(impatch, fprofile, z);
            else
                goodpatch = true;
            end
            
            % Assign patch to new stack
            if goodpatch
                J(:,:,i) = impatch;
            else
                continue
            end
        end

        % Warning and return when too many resamples performed for a single
        % patch due to low signal present in patch
        if resamplecount >= resamplelimit
            wmsg = ['Please consider increasing the window size (fov)' ...
                    newline...
                    'Small patch sizes may tend to have a lower ' ...
                    'average SNR or will be more likely to resample ' ...
                    'areas with low SNR due to sparse signal.' ...
                    newline...
                    'Only ' num2str(i) ' out of ' num2str(numpatches)...
                    ' patches were successfully sampled'];
            warning(wmsg)
            return
        end
    end
end


function isgood = CheckSNR(I, fprofile, z)
% CHECKSNR approximately checks that the signal is sufficiently high in the
% image slice to consider it a valid or "good" representative sample.
    actual   = mean(I(:));
    expected = fprofile(z);
    isgood   = actual >= expected;
end