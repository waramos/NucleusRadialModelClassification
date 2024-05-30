function C = FindBlobs(I, s1, s2, thresh, r)

    % DoG
    I = SimpleDoG(I, s1, s2);

    % Threshold on DoG Image
    I = ThresholdCleanly(I, thresh, r);

    % Get watershedded mask
    I = SimpleWaterShed(I);

    % Blobs are centroids of remaining regions
    C = regionprops(I, 'Centroid');
    C = cat(1, C.Centroid);
end