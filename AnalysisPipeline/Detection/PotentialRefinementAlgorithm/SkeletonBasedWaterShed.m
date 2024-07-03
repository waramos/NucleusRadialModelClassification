function P = SkeletonBasedWaterShed(Mask)
% SKELETONBASEDWATERSHED will approximate a skeleton of a binary image and
% then proceeed to pull points that are approximate centroids of 

    % Pseudo-branching requirement (?) - i.e. how many pixels need to be on
    % in the local nieghborhood and size of neighborhood (w x w window)
    w  = 3;
    h  = ones(w);
    nv = numel(h)-w;

    % Skeletonizes mask
    BW = bwdist(~Mask);
    MF = medfilt2(BW);
    S  = BW>MF;

    % Preserve the approx. centroids and discard other skel regions
    od = ordfilt2(S, nv, h);
    P  = regionprops(od, 'Centroid');
    P  = cat(1, P.Centroid);
end