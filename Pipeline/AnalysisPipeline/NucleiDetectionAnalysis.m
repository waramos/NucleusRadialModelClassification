% An analysis of the spread of points identified within each cluster as a
% result of using a given radius, r, in DBSCAN...
% Rather than using r = r_NN, r = dz_px or even r = sqrt(dz^2+(r_NN/2)^2) 
% might be better in terms of clustering only a single cell's detections.
% This might increase the separation. 
addpath('D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\IlluminationProfileCorrected\Results\')
load NuclearDetections.mat


% iso scaling?
isoscaling = true;

for i = 1:3
    C       = Nuclei(i).Centroids;
    OP      = Nuclei(i).OriginalPoints;
    labels  = Nuclei(i).Cluster;
    nlabels = numel(unique(labels));
    % Loop over clusters
    for c = 1:nlabels
        cell_centroid = C(c, :);

        % index PC
        c_idx = labels == c;
        [c_idx, ~] = find(c_idx);
        c_pc  = OP(c_idx, :);
        c_pc  = c_pc - cell_centroid;
        OP(c_idx, :) = c_pc;
    end

    if isoscaling
        OP = OP*dx;
    end

    OP_c{i} = OP;
end


%% Visualizations

figure; 
plot3(OP_c{1}(:,1), OP_c{1}(:,2), OP_c{1}(:,3), '.r')
xlabel('X'); 
ylabel('Y'); 
zlabel('Z')
axis image

[coeff, score, latent] = pca(OP_c{1});


[theta, rho, z] = cart2pol(OP_c{1}(:,1), OP_c{1}(:,2), OP_c{1}(:,3));

cv_AngleR = cov([theta, rho, z]);
cv_XYZ = cov(theta, rho);

% figure; imagesc(cv_AngleR)

