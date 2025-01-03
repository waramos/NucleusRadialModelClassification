function P = SetPointsCoordinateSys(Tbl, xypxsz, zspacing, ispixelunits)
% SETPOINTSCOORDINATESYS will extract 3D point cloud data from a table and
% then transform it to be an isotropic pointcloud in either pixel units or
% microns. The units used for the isotropic coordinate transformation are
% determined by a logical.
%
% INPUTS:
% T - data table of point detections across 2D slices in a volume
% xypxsz - x / y pixel sizing across image in microns
% zspacing - image volume spacing between planes, in microns
% ispixelunits - false if user wants micron units, true if pixel units
% desired
%
% OUPUTS:
% P - pointcloud transformed to be isotropic and in the desired coordinate
% system


    % Logical for setting the point cloud in pixel coordinates or physical
    % coordinate space
    if nargin < 4 || isempty(ispixelunits)
        ispixelunits = true;
    end

    hasZslices = any(contains(Tbl.Properties.VariableNames, 'Z'));
    hasTpoints = any(contains(Tbl.Properties.VariableNames, 'T'));

    % Pointcloud data extracted from the table
    X = Tbl.X;
    Y = Tbl.Y;
    
    % Has Z slices
    if hasZslices
        Z = Tbl.Z;
    end

    % Has timepoints/frames
    if hasTpoints
        T = Tbl.T;
    end

    % Scaling z slices isotropically - units: in pixels
    if hasZslices
        zfactor = (zspacing/xypxsz);
        Z       = Z*zfactor;
    end
    
    % Microscope/real coordinate space - units: in microns
    if ~ispixelunits
        X = X*xypxsz;
        Y = Y*xypxsz;
        if hasZslices
            Z = Z*xypxsz;
        end
    end

    % Arranges coordinates in columns so each column is a coord dimension
    P = [X Y];
    
    if hasZslices
        P = [P Z];
    end

    if hasTpoints
        P = [P T];
    end
end