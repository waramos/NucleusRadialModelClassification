function [idz, d_df, df, fprofile] = EstimateFocalPlaneIndices(I)
% ESTIMATEFOCALPLANEINDICES will estimate which z slices from an image
% volume are in focus. This is done by first computing the average signal
% across each plane. Then, the gradient of the average fluorescence with
% respect to z is computed. This is smoothed and the point at which the
% greatest gradient increase occurs is approximated to be where content
% first most noticeably comes into focus.
    
%% Fluorescence change
    % Average signal profile across the stack
    fprofile = mean(I, [1 2]);
    fprofile = squeeze(fprofile);
    
    % dF/dz
    df = conv(fprofile, [1 0 -1]/2, 'same');
    df = smooth(df(2:end-1));

    % Approx inflection point where things begin to be in focus
    d_df     = conv(df, [1 0 -1]/2, 'same');
    [~, idz] = max(d_df(2:end-1));

    % Adjustment of index per clipping of gradients
    idz = idz + 2;
end