function Mask = HardThreshold(I, threshold)
% HARDTHRESHOLD will median filter an image and then apply a scalar value
% as a threshold value. This function was adapted from the "Simple
% Intensity Threshold" method defined in the segmentation GUI to consider
% higher dimensions (3-4 for now) to operate on quickly.


    % Looping through the third and fourth dimensions of the data
    [D3, D4] = size(I, [3 4]);

    for j = 1:D4
        for i = 1:D3
            I(:, :, i, j) = medfilt2(I(:,:,i, j), [3 3], 'symmetric');
        end
    end

    disp('Median Filtering Done')

    % Simple threshold of the median filtered image
    Mask = I>threshold;
end