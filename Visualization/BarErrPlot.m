function BarErrPlot(dmat, errmat)
% BARERRPLOT will plot a bar graph with errors given that the user has
% passed in a data matrix, dmat, and its corresponding error matrix,
% errmat.

    figure
    bp = bar(dmat');
    hold on

    % X coordinates of the bars
    xc = vertcat(bp.XEndPoints);

    % Error plot atop the bars
    errorbar(xc', dmat', errmat', 'k','linestyle','none');

end