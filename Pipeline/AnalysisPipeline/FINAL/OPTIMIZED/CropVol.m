function Vol = CropVol(Vol, cropregion)
% CROPVOL will crop a volume based on the region specified by the
% cropregion matrix. First row of matrix is starting index and second row
% is the ending index.

    % Get crop coordinates
    x   = cropregion(1,1):cropregion(2,1);
    y   = cropregion(1,2):cropregion(2,2);
    z   = cropregion(1,3):cropregion(2,3);
    Vol = Vol(y, x, z);
end