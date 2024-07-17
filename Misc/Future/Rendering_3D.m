ImVol = I(1:1024, 1:1024, :);
f = uifigure;
VV = viewer3d(f);
VS = volshow(ImVol, 'Parent', VV);
VV.RenderingQuality = 'high';
dx = 6.5/25; dy = dx; dz = .5;
v    = ones(4,1);
v(1) = dx;
v(2) = dy;
v(3) = dz;
v    = diag(v);
VS.Transformation.A = v;
C = Nuclei(end).Centroids;
C(:,3) = C(:,3)*rsfactor;
idx_c = (C(:,1) <= 1024) & (C(:,2)<=1024);
rows = find(idx_c);
C2 = C(rows, :);
BinVol = false(1024, 1024, 80);
C3 = round(C2);
C3 = max(C3, 1);
C3 = min(C3, 1024);
for i = 1:size(C2);
y = C3(i,1);
x = C3(i,2);
z = C3(i,3);
BinVol(x, y, z) = true;
end
BW = bwdist(BinVol)<6;
VS.OverlayData  = BW;
VS.GradientOpacityValue
VS.GradientOpacityValue = 0.1
VS.GradientOpacityValue = 1
VS.GradientOpacityValue = 0.1
VS.RenderingStyle
VS.RenderingStyle = 'GradientOpacity';
VS.GradientOpacityValue = 0.2;
clc
VS.GradientOpacityValue = 0.5;
VS.GradientOpacityValue = 0.01;
VS.GradientOpacityValue = 0.05;
VS.AlphaData
VS.Alphamap
VS.Alphamap = linspace(0.25, 0.75, 256);
VS.Alphamap = linspace(0.25, 0.3, 256);
VS.Alphamap = linspace(0.55, 1, 256);
VS.Alphamap = linspace(0, .5, 256);
VS.Alphamap = linspace(0, .75, 256);

cv = linspace(0, 1, 256);
cmap = [cv' zeros(256, 1) cv'];

VS.Colormap = cmap;
VV.BackgroundColor = [1 1 1];
VV.BackgroundGradient = false;

VS.Alphamap = linspace(0, .1, 256);
VS.GradientOpacityValue = 0.04;
VS.Alphamap = linspace(0, .25, 256);