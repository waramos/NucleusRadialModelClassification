function Mask = ThresholdCleanly(I, thresh, r)
    Mask = I>thresh;

    % Disk kernel
    h = -r:r;
    h = h.^2;
    h = h+h';
    h = sqrt(h)<=r;

    % Morph open
    Mask = imopen(Mask, h);

    % Morph close
    Mask = imclose(Mask, h);
end