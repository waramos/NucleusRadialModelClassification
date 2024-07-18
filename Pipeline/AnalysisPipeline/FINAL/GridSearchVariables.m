function params = GridSearchVariables(ipvalues, stepsz, numsteps)
% GRIDSEARCHVARIABLES will produce an MxMx2 grid of values to search across
% given initial parameter values, a step size, and number of steps. The
% first index in the third dimension will leave values across rows
% (dimension 1) the same but change across columns (dimension 2). The
% second index in the third dimension will instead do the converse with
% respect to its values.
%
% INPUTS:
%
% ipvalues (1x2 row vector) - initial parameter values
%
% stepsz (scalar) - step size between values
%
% numsteps (scalar) - number of steps / values to try, including the
% initial values
%
%
% OUTPUTS:
%
% params (MxMx2 matrix) - grid of values to try.
%
%
% William A Ramos, Kumar Lab @MBL July 2024


    % Hard coded for just two variables for now - can consider generalizing
    % later
    nparams     = 2;
    
    % Final value after n steps
    finalvalues = ipvalues + stepsz*(numsteps-1);
    pvalues     = zeros(numsteps, 2);
    for i = 1:nparams
        pvalues(:,i) = (ipvalues(i):stepsz:finalvalues(i))';
    end

    % Create the mesh of values to use for the parameter experiment
    [X, Y] = meshgrid(pvalues(:,1), pvalues(:,2));
    params = cat(3, X, Y);
end