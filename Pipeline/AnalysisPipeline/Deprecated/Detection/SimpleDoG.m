function J = SimpleDoG(I, sigma1, sigma2)
% SIMPLEDOG is a simple difference of gaussians.
    I1 = imgaussfilt(I, sigma1);
    I2 = imgaussfilt(I, sigma2);
    J  = I1-I2;
end