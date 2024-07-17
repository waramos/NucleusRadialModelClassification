function J = IsoResampleImageStack(I, dx, dz, memoryeff)
% RESAMPLEIMAGESTACK will efficiently resample an image stack. By default,
% this function is optimized to use minimal memory but if the user has
% sufficient memory and prefers computational efficiency instead, the
% function will permute the array to interpolate along z and then permute
% the dimensions back.
%
% William A. Ramos, Kumar Lab @ MBL - July 2024

    if nargin < 4 || isempty(memoryeff)
        memoryeff = true;
    end
    
    % Getting dimensions to resample
    [M, N, Z] = size(I, [1 2 3]);

    % Figuring out number of new z slices
    dz_dx = dz/dx;
    Z2    = floor(Z*dz_dx);

    % Original sample points
    zs = 1:dz_dx:Z2;
    zq = 1:Z2;

    % Resampling on a row by row basis
    ogclass = class(I);
    % Memory efficient
    I       = single(I);

    % The function will branch in terms of how it interpolates depending on
    % whether the user wants to be memory efficient
    if memoryeff
        % Init array to assign to
        J       = zeros(M, N, Z2, 'single');
        for r = 1:M
            % Pull out the first row
            Im = I(r, :, :);
            
            % Faster when only operating on rows
            Im = squeeze(Im);
            Im = interp1(zs, Im', zq);
            Im = Im';
            Im = reshape(Im, 1, N, Z2);
            
            % Assignment on a row by row basis
            J(r,:,:) = Im;

            % Reporting on progress
            pct      = (r/M)*100;
            msg      = ['Percent resampled: ' num2str(pct, '%.2f') '%'];
            disp(msg)
        end
    else
        % If memory is not a concern, then the user can easily interpolate
        % the entire stack 
        I = permute(I, [3 1 2]);
        J = interp1(zs, I, zq);
        J = permute(J, [2 3 1]);
    end

    % Cast at end for efficiency
    J = cast(J, ogclass);
end