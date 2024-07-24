function h = MakeDiskKernel(r)
% MAKEDISKKERNEL will make a 2D disk kernel based off of a given radius, r.
    h = -r:r;
    h = h.^2;
    h = h+h';
    h = sqrt(h)<=r;
end