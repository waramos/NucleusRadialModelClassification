function [C, I] = BlobsWithSkelWaterShed(I, s1, s2, thresh, r)
% FINDBLOBS will return both the reduced cluster points and image mask
% based on a difference of gaussian blob detection approach.

    % DoG
    I = SimpleDoG(I, s1, s2);

    % Threshold on DoG Image with morphological filtering and flood fill
    I = ThresholdCleanly(I, thresh, r);

    % Get watershedded mask
    C = SkeletonBasedWaterShed(I);
end