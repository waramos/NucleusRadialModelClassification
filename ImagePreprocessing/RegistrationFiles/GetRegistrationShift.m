function [xshift, yshift] = GetRegistrationShift(I, J)
% GETREGISTRATIONSHIFT will register two images using phase correlation and
% get subpixel registration by performing a gaussian fit on an 11x11 pixel
% window near the integer pixel peak of the power spectrum. This function
% makes the assumption that there is no z shift in the data.

% Phase correlation estimation of shifts between slices
    numzslices = size(I, 3);
    Q          = zeros(size(I, [1 2]), 'single');
    for z = 1:numzslices

        % Pulling out slices from the 
        im1 = I(:, :, z);
        im2 = J(:, :, z);
    
        % Denoising
        im1 = medfilt2(im1);
        im2 = medfilt2(im2);

        % Z sum of the phase correlation to
        q   = FourierShift2D(im1, im2);
        Q   = Q+q;
    end
    
    [xs, ys]   = SubPixelFit(Q);


        function [x, y] = SubPixelFit(obj, Q)
            % SUBPIXELFIT will perform a gaussian fit on a window of the
            % summed power spectrum near its global maxima to find a sub
            % pixel registration.

            % Fitting window in pixels (+/-)
            wdw = 5;
            sz = size(Q);

            [~,midx] = max(Q,[],'all');
            [yi,xi]  = ind2sub(size(Q),midx);
            [xl,yl]  = BoundCheck([xi-wdw xi+wdw],[yi-wdw yi+wdw],sz);
            
            grab               = Q(yl(1):yl(2),xl(1):xl(2));
            [fparams,sparams]  = obj.FitGauss(grab,[xl(1) yl(1)],'lsq');
            
            if ~isempty(fparams) % Find x and y shifts in original coordinate system
                x = fparams(1,2)-sz(2)/2-1+mod(sz(2),2)*0.5;
                y = fparams(1,3)-sz(1)/2-1+mod(sz(1),2)*0.5;
            else
                x = sparams(1,2)-sz(2)/2-1+mod(sz(2),2)*0.5;
                y = sparams(1,3)-sz(1)/2-1+mod(sz(1),2)*0.5;
                disp('Fit unsuccessful: Reverted to centroid.')
            end

            % Bound checks for grab box
            function [xl,yl] = BoundCheck(xl,yl,sz)
                if xl(1)<1;     xl = xl+abs(min(xl))+1; end
                if yl(1)<1;     yl = yl+abs(min(yl))+1; end
                if xl(2)>sz(2); xl = xl+(sz(2)-xl(2));  end
                if yl(2)>sz(1); yl = yl+(sz(1)-yl(2));  end                
            end
        end