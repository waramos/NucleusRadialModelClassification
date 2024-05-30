function W = SimpleWaterShed(Mask)
    BW       = bwdist(~Mask);
    D        = -BW;
    W        = watershed(D);
    W(~Mask) = 0;
    W        = W>0;
end