
% Getting all of the Expression Result files
imds   = imageDatastore('D:\MichaelPalmer_Gillis\2024_04_26 - ATDC5 Cells\IlluminationProfileCorrected\Results\SensitivityAnalysis', 'FileExtensions', '.mat');
fNames = imds.Files;
nfiles = numel(fNames);
ncond  = sqrt(nfiles);

% Creating the two figures to show all results in a grid like layout
fig1          = figure;
fig1.Units    = "normalized";
fig1.Position = [0.05 0.05 0.45 0.85];
fig2          = figure;
fig2.Units    = "normalized";
fig2.Position = [0.5 0.05 0.45 0.85];

fNames2 = cat(1, fNames(6:end), fNames(1:5));

LeftMost   = 1:5:21;
BottomMost = 21:25;

% Producing the grid layout of plots for direct comparison
for i = 1:nfiles

    % Pull out file and load results
    fid = fNames2{i};
    load(fid, 'ExpressionResults')
    [~, fn, ~] = fileparts(fid);

    % Create new axes
    ax1 = subplot(ncond, ncond, i, 'Parent', fig1);
    ax2 = subplot(ncond, ncond, i, 'Parent', fig2);

    % Plot directly atop the new axes
    % idx1 = regexp(fn, 'P1');
    % idx2 = regexp(fn, '_results');
    % tmod = fn(idx1:idx2-1);

    if any(i == LeftMost)
        ylab = true;
    else
        ylab = false;
    end

    if any(i == BottomMost)
        xlab = true;
    else
        xlab = false;
    end

    tmod = '';
    [newStats, ALLRESULTS] = GetStatsAndPlots(ExpressionResults, tmod, ax1, ax2, xlab, ylab);
    Stats{i,1} = newStats;
    Stats{i,2} = ALLRESULTS;

    % updates the graphics
    drawnow
end