function J = SimpleDoG(I, sigma1, sigma2)
    I1 = imgaussfilt(I, sigma1);
    I2 = imgaussfilt(I, sigma2);
    J  = I1-I2;
end