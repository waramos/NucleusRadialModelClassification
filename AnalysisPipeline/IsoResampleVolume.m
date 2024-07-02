function V2 = IsoResampleVolume(V, dxdy, dz)
% ISORESAMPLEVOLUME will resample a volume to be isotropic if it is not
% already isotropic.

    % Throw error for isotropic cases
    if dxdy == dz
        V2 = V;
        warning(['Resolution is already isotropic per specified'...
              ' xy pixel size and z spacing. No computation performed.'])
        return
    end

    % Scaling z pixels accordingly
    dz_dx = dz/dxdy; 

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

    % Interpolation method
    if islogical(V)
        rsmethod = 'nearest';

    elseif isa(V, 'uint16') || isa(V, 'uint8')
        rsmethod = 'linear';

    else
        rsmethod = 'makima';

    end

    % Resampled image
    V2 = interp3(Xs, Ys, Zs, V, Xq, Yq, Zq, rsmethod);
end
