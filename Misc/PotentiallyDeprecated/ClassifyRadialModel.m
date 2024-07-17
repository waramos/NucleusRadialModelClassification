function [O, BW, CCData, localizationerr] = ClassifyRadialModel(P, r, V, voxelsz, p_iso)
% CLASSIFYRADIALMODEL looks at a pointcloud, P, to see if there is overlap
% between spheres at P with radius, r, and a binary 3D volume, V. We assume
% that P and V have the same voxel resolution, i.e. P originates from
% detections in an image stack with the same dimensionality and voxel
% sizing as V, but isotropic scaling can be done if needed.

    if nargin < 5 || isempty(p_iso)
        p_iso = true;
    end

    % Extracting voxel information in case of anisotropy
    if nargin < 4 || isempty(voxelsz)
        voxelsz = [1 1 1];
    end

    % Voxel resolution in x, y, and z dimensions
    dx = voxelsz(1);
    dy = voxelsz(2);
    dz = voxelsz(3);
    
    if dx ~= dy
        % Expects isotropic resolution in xy planes
        error('X and Y are not isotropic...')
    else
        % conversion factor: z in # of xy pixels
        if ~p_iso
            dz_dx = dz/dx; 
        else
            dz_dx = dz;
        end
    end

    % Dimension information from channel of interest
    [M, N, Z] = size(V, [1 2 3]);
    Zpx       = Z*dz_dx;          % last z index * conversion factor

    % Sampling vectors
    x     = 1:N;
    y     = 1:M;
    z_px  = 1:dz_dx:Zpx; % sampled:  xy pixel res scaled z slices 
    z_px2 = 1:Zpx;       % querying: integer z slicing for new sampling
    Zpx   = z_px2(end);

    % Interpolation Mesh for ~isotropic resolution
    [Xs, Ys, Zs] = meshgrid(x, y, z_px);  % sample points
    [Xq, Yq, Zq] = meshgrid(x, y, z_px2); % query points for new stacks

    % Interpolated channel of interest
    newV = interp3(Xs, Ys, Zs, V, Xq, Yq, Zq, 'nearest'); %

    % Interpolated reference channel from point clouds
    Mask = false(size(newV));

    % Scaling by voxel resolution if not isotropic
    if p_iso
        P2 = P;       % if already isotropic, leave as is
    else
        P2 = P*dz_dx; % xy pixel sampling in z
    end

    % Bound checking
    P2        = max(1, P2);
    P2(:,1)   = min(P2(:,1), N);
    P2(:,2)   = min(P2(:,2), M);
    P2(:,3)   = min(P2(:,3), Zpx);
    P2        = round(P2);    

    % Isotropic Dilation to create Spheres in Volume
    idx       = sub2ind([M N Zpx], P2(:,2), P2(:,1), P2(:,3));
    Mask(idx) = 1;
    BW        = bwdist(Mask, 'euclidean');
    BW        = BW < r;

    % Converting the volume to voxel points
    idx          = find(V);
    [x, y, z_px] = ind2sub([M N Zpx], idx);

    % Overlap
    O = bwselect3(BW, x, y, z_px, 26);

    % MSE from subpixel
    SE                    = (P - P2).^2;
    localizationerr.MSE   = mean(SE, 1);
    localizationerr.SE    = SE;          % gives info on most/least off

    % Connected components
    CCData.Overlapping = bwconncomp(O);
    CCData.Original    = bwconncomp(BW);
end