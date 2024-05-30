function V = QuickBeamCorrection(V)

    % Correcting the beam artifact at edges
    Stackmean1  = mean(V, [1 2]);
    ogclass     = class(V);
    Lineprofile = mean(V, [1 3]);
    Lineprofile = squeeze(Lineprofile);
    Profile     = max(Lineprofile(:))./Lineprofile;
    V           = single(V);
    V           = V.*Profile;
    Stackmean2  = mean(V, [1 2]);

    % Scaling according to magnitude difference between normalized image
    % and the orignal image volume
    Stackmean1  = single(Stackmean1);
    rsfactor    = Stackmean1./Stackmean2;
    V           = V.*max(rsfactor);
    V           = cast(V, ogclass);
end