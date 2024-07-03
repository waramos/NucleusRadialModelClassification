function mdata = GetH5ImageParams(fid)
% GETH5IMAGEPARAMS will get the imaging parameters used to acquire the
% data. This assumes that the input is a filepath to an .h5 image which
% shares filename with a .tif file.

    % Voxel information
    pixelsize = h5read(fid, '/pixelsize');
    dz        = h5read(fid, '/DzStage');

    % Organizing into data struct
    mdata.pixelsize  = pixelsize;
    mdata.zstep      = dz;
    mdata.voxelsize  = [pixelsize pixelsize dz]; % Assumes 2D, xy isotropic
    mdata.voxelunits = 'microns';                % Default assumption
end