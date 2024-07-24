function Mask = RefineMask(Mask, Radius)
    if Radius > 0
        % Closes holes and improves boundary of binary image
        Radius = round(Radius);
        
        % Circular mask kernel
        d    = 2*Radius+1;
        cx   = 1:d;
        cy   = cx';
        cx   = (cx-(Radius+1)).^2;
        cy   = (cy-(Radius+1)).^2;
        h    = (cx + cy) <= ((Radius)^2);

        % Morphological clean up operations
        Mask = imopen(Mask, h);
        Mask = imdilate(Mask,[0 1 0;1 1 1;0 1 0]);

        disp('Morphological Filtering Done')
    end
end