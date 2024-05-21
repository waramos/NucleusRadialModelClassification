function J = BootstrapImageVolume(I, mrows, ncols, numpatches)
% BOOSTRAPDATA will bootstrap a 3D image stack with some field of view
% based on number of specified rows and columns. It will produce samples as
% patches with the requested field of view and as many samples as are
% requested.

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

    % Init shuffled image array
    ogclass = class(I);
    J       = zeros(mrows, ncols, numpatches, ogclass);
    for i = 1:numpatches
        % Start and end indices
        y1 = idy(i);
        x1 = idx(i);
        y2 = y1 + mrows - 1;
        x2 = x1 + ncols - 1;

        % Patch indices
        y  = y1:y2;
        x  = x1:x2;
        z  = idz(i);

        % Assign patch to new stack
        J(:,:,i) = I(y, x, z);
    end
end