function [Vols, radii] = GetRadialShells(r, dr)

    % Range of radial distances to consider
    if numel(r) == 2
        r1 = r(1);
        r2 = r(2);
    else
        r1 = 0;
        r2 = r;
    end

    % Spherical Volume equation
    Vol   = @(r) 4*pi*r.^3;

    % Shells' inner radii and their respective volumes
    radii = r1:dr:r2;
    Vols  = Vol(radii);
    Vols  = Vols(2:end)-Vols(1:end-1);
end