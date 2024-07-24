function [y, x, z] = GetLocalCellNHood(pts, r, sz)
% GETLOCALCELLNHOOD constructs a local neighborhood based on a pointcloud
% of detections and returns the indices corresponding to 3 dimensions.
    [x1, x2] = bounds(pts(:,1));
    [y1, y2] = bounds(pts(:,2));
    [z1, z2] = bounds(pts(:,3));

    % X Y Z coordinates for the local neighborhood/window to grab
    y = (y1-r):(y2+r);
    x = (x1-r):(x2+r);
    z = (z1-r):(z2+r);
    y = round(y);
    x = round(x);
    z = round(z);

    % Bound check
    y = min(y, sz(1));
    x = min(x, sz(2));
    z = min(z, sz(3));
    y = max(y, 1);
    x = max(x, 1);
    z = max(z, 1);
end