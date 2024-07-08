function Rstats = GetNuclearDist(P)
% GETNUCLEARRADIUS will get the approximate radius as well as other
% information from the point cloud representing nuclei detections.

    % Number of z slices and timepoints
    numz = max(P(:,3));
    numt = max(P(:,4));

    % Distance to nearest neighbor
    Rstats = struct('Dist2NN', []);

    % Init tracking var
    count     = 0;
    numimages = numt*numz;

    for t = 1:numt
        % Indices for the given timepoint
        idt = P(:,4) == t;
        for z = 1:numz
            % Increase count
            count = count + 1;

            % Indices for the given z slice
            idz   = P(:,3) == z;

            % Time and z indices combined
            idx   = idt & idz;
            r     = find(idx);
            Cn    = P(r, 1:2);

            % No points for a given 2D image goes to next iteration
            if isempty(Cn) || size(Cn, 1) < 2
                continue
            end

            % Finding centroids' nearest neighbor
            [~, D] = knnsearch(Cn, Cn, 'K', 2);
            C_l2   = D(:,2);

            % Assigning the nearest-neighbor distance
            Rstats(z, t).Dist2NN = C_l2;

            % CL reporting
            pct   = count/numimages;
            pct   = pct*100;
            msg   = ['Percent computed: ' num2str(pct, '%.2f')];
            disp(msg)
        end
    end
end