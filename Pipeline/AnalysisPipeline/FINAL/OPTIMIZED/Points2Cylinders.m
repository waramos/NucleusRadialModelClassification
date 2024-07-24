function Nuclei = Points2Cylinders(pts, r, sz)

    % Round to nearest pixel - new var enables comparison
    PC_pix = round(pts);
    x      = PC_pix(:,1);
    y      = PC_pix(:,2);
    z      = PC_pix(:,3);

    % Convert point cloud to linear indices corresponding to pixels
    idx = sub2ind(sz, y, x, z);

    % Init nuclei array and set the true values
    Nuclei      = false(sz(1), sz(2), sz(3));
    Nuclei(idx) = true;

    % Convolve Nuclei stack with disk kernal for cylinder model
    h      = MakeDiskKernel(r);
    Nuclei = imdilate(Nuclei, h);
end