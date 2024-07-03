function V = IsoResampleVolume(V, dxdy, dz)
% ISORESAMPLEVOLUME will resample a volume to be isotropic if it is not
% already isotropic.

    % Throw error for isotropic cases
    if dxdy == dz
        warning(['Resolution is already isotropic per specified'...
              ' xy pixel size and z spacing. No computation performed.'])
        return
    end

    % Scaling z pixels accordingly
    dz_dx = dz/dxdy; 

    % Dimension information from channel of interest
    [M, N, Z] = size(V, [1 2 3]);
    Z2        = ceil(Z*dz_dx);          % last z index * conversion factor

    % Sampling vectors
    x     = 1:N;
    y     = 1:M;
    z_px  = 1:dz_dx:Z2; % sampled:  xy pixel res scaled z slices 
    z_px2 = 1:Z2;       % querying: integer z slicing for new sampling

    % Interpolation Mesh for ~isotropic resolution
    [Xs, Ys, Zs] = meshgrid(x, y, z_px);  % sample points
    Xs           = single(Xs);
    Ys           = single(Ys);
    Zs           = single(Zs);

    % Query points for new stacks
    [Xq, Yq, Zq] = meshgrid(x, y, z_px2); 
    Xq           = single(Xq);
    Yq           = single(Yq);
    Zq           = single(Zq);

    % Interpolation method
    if islogical(V)
        rsmethod = 'nearest';
    else
        rsmethod = 'linear';
    end

    % Resampled image
    V  = single(V);
    if numel(y) > 1
        % When there is some y information to consider
        V = interp3(Xs, Ys, Zs, V, Xq, Yq, Zq, rsmethod);
    else
        % Faster when only operating on rows
        V = squeeze(V);
        V = interp1(z_px, V', z_px2);
        V = V';
        V = reshape(V, M, N, Z2);
    end
end
