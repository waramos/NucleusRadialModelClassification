function [Majmom, Minmom] = MomentsOfInertia(x, y, datadir)
    
    if nargin < 3 || isempty(datadir)
        datadir = true;
    end

    % "Centroid" of data
    c      = [mean(x) mean(y)];
    
    % Moments of inertia
    C      = zeros(2);
    dx     = x-c(1);
    dy     = y-c(2);
    C(1,1) = mean(dx.^2);
    C(2,1) = mean(dx.*dy);
    C(1,2) = C(2,1);
    C(2,2) = mean(dy.^2);
    [V,D]  = eigs(C);
    D      = 2*sqrt(diag(D));
    
    % Major and minor moments / eigen vecs as coordinates
    AxesCoords = [-V(:,1)'; V(:,2)'];
    Majmom     = [c(:)';c(:)'+D(1)*AxesCoords(1,:)];        
    Minmom     = [c(:)';c(:)'+D(2)*AxesCoords(2,:)];

    

    % Axes go in both directions rather than solely in eig dir
    % if ~datadir
        % Span across both dimensions
        Majdx       = diff(Majmom(:,1));
        Majdy       = diff(Majmom(:,2));
        Majmom(1,1) = Majmom(1,1) - Majdx;
        Majmom(1,2) = Majmom(1,2) - Majdy;
        Mindx       = diff(Minmom(:,1));
        Mindy       = diff(Minmom(:,2));
        Minmom(1,1) = Minmom(1,1) - Mindx;
        Minmom(1,2) = Minmom(1,2) - Mindy;
    % end
end