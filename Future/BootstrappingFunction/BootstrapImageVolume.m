function [J, Indices] = BootstrapImageVolume(I, patchsize, numpatches, maxproj, zs, fprofile, resamplelimit)
% BOOSTRAPIMAGEVOLUME bootstraps a 3D image stack with some field of view
% based on number of specified rows and columns. It will produce samples as
% patches with the requested field of view and as many samples as are
% requested.
% 
% INPUTS:
% I (3D array) - image stack to perform random patch grabbing on
%
% patchsize (3 element vector) - specifies how many rows, columns, and z 
% slices to grab for a patch
%
% numpatches (scalar) - specifies how many 3D patches/slabs the user would
% like to sample 
%
% maxproj (logical) - specifies whether to max project 3D patches along the
% third dimension. Will do so if true. True by default.
%
% zs (2 element vector) - starting and ending index along z in case there 
% is information that is out of focus outside the window specified by zs(1)
% and zs(2).
%
% fprofile (vector) - column vector containing the mean fluorescence /
% pixel intensity of 2D slices across z
%
%
% OUTPUTS:
% J (3D array) - new image stack that will be MxNxZ where M and N are the
% first two elements of the patch size vector and Z will either be equal to
% numpatches (if maxproj was true or patchsize(3) is 1) or
% numpatches*patchsize(3) if maxproj is false.
%
% Indices (data struct) - contains the indices that were randomly selected,
% patch size information, if maximum projection was requested, number of 
% requested patches, and the z window used.
%
% William A. Ramos, Kumar Lab @ MBL July 2024

    % Default patch size is 3D slab
    if nargin < 2
        patchsize = [64 64 3];
    end

    % Number of patches to randomly grab from the image stack
    if nargin < 3 || isempty(numpatches)
        numpatches = 10;
    end

    % Performs max projection of 3D by default
    if nargin < 4|| isempty(maxproj)
        maxproj = true;
    end

    % Starting z slice where information begins to be in focus
    if nargin < 5 || isempty(zs)
        zs = 1;
    end

    % Fluorescence profile across the stack for consideration of sampling
    if nargin < 6 || isempty(fprofile)
        fprofile = [];
    end

    % When fprofile is known and only images with in focus content are
    % desired, there is a limit to resampling patches to eventually acquire
    % a patch that is in focus
    if nargin < 7 || isempty(resamplelimit) 
        resamplelimit = 1000;
    end

    % Patch size info
    mrows   = patchsize(1);
    ncols   = patchsize(2);
    zslices = patchsize(3);

    % Determines whether user is wanting a 3D patch
    if zslices > 1
        patch3D = true;
    end

    % Only use planes known to be in focus
    z1 = zs(1);
    if numel(zs) == 2
        z2 = zs(2);
    else
        z2 = size(I, 3);
    end

    % Z windowing/cropping
    I = I(:,:,z1:z2);

    % Minimum window sizing
    mn  = 64; % min rows / cols
    mnz = 1;  % min z slices - ensures this works on 2D data

    % Get image stack bounds
    [M, N, Z] = size(I, [1 2 3]);

    % 2D consideration
    if Z == 1
        zslices = 1;
    end

    % When % of dimension requested - ensure integer valued FOV sizing
    if mrows < 1
        mrows = round(mrows*M);
    end
    if ncols < 1
        ncols = round(ncols*N);
    end
    if zslices < 1
        zslices = round(zslices*Z);
    end

    % Bound checking sets upper limit on randomization of starting index
    mrows   = max(mrows, mn);
    ncols   = max(ncols, mn);
    zslices = max(zslices, mnz);
    newM    = M - mrows;
    newN    = N - ncols;
    newZ    = Z - zslices;
    newZ    = max(newZ - zslices, 1);

    % Sampling indices - i.e. starting index of random patch
    idy = randi(newM, [numpatches 1]);
    idx = randi(newN, [numpatches 1]);
    idz = randi(newZ, [numpatches 1]);

    % Indices info for user to know where patches came from
    Indices.XCols      = idx;
    Indices.YCols      = idy;
    Indices.ZSlices    = idz;
    Indices.PatchSize  = patchsize;
    Indices.MaxProject = maxproj;
    Indices.NSamples   = numpatches;
    Indices.ZWindow    = zs;

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
            z1 = idz(i);
            y2 = y1 + mrows - 1;
            x2 = x1 + ncols - 1;
            z2 = z1 + zslices - 1;
    
            % Patch indices
            y  = y1:y2;
            x  = x1:x2;
            z  = z1:z2;
    
            % Grabbing image patch;
            impatch = I(y, x, z);
            % Maximum projection option
            if maxproj
                impatch = max(impatch, [], 3);
            end

            % Compare mean intensity of current patch and the z index
            if ~isempty(fprofile)
                % Adjust for dropping z indices
                z         = min(z-z1);
                goodpatch = CheckSNR(impatch, fprofile, z);
            else
                goodpatch = true;
            end
            
            % Assign patch to new stack
            if goodpatch
                % Optional 3D or 2D patches
                if patch3D && ~maxproj
                    nz       = size(impatch, 3);
                    z1       = (i-1)*nz + 1;
                    z2       = z1 + nz;
                    z        = z1:z2;
                else
                    z = i;
                end
                J(:,:,z) = impatch;
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

    % Max projects the patch
    if size(I, 3) > 1
        I = max(I, [], 3);
    end

    actual   = mean(I(:));
    expected = fprofile(z);
    isgood   = actual >= expected;
end