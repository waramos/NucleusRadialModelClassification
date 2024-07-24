function PC = AdjustPointCloud(PC, sz, rsfactor)
% ADJUSTPOINTCLOUD ensures the point cloud is rounded to the nearest pixel,
% then corrected by the rescale factor, and lastly bound checks the final
% point cloud to make sure it can index the image accurately.
    
    % Only rescale if the value is not 1
    if rsfactor ~= 1
        PC(:,3) = PC(:,3) * rsfactor;
    end

    % Converting size vector to x y z ordering
    sz = [sz(2) sz(1) sz(3)];
    
    % Bound check the point cloud
    PC = round(PC);
    PC = max(1, PC);
    PC = min(PC, sz);
end