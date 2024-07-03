function W = SimpleWaterShed(Mask)
% SIMPLEWATERSHED applies a conventional watershed to masks to try to
% separate out blobs into different connected components.
    BW       = bwdist(~Mask);
    D        = -BW;
    W        = watershed(D);
    W(~Mask) = 0;
    W        = W>0;
end