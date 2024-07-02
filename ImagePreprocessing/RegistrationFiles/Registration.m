% Image registration class
% [MP]
classdef Registration < handle
    properties (Access=public)
        verbose   (1,1) logical = 0  % Verbose mode on/off
        averaging (1,1) logical = 1  % Average through unused dimensions?
        pixelmode char          = 'pixel'
    end
    
    methods
        function set.pixelmode(obj,val)
            pm = lower(val);
            if strcmp(pm,'pixel') || strcmp(pm,'subpixel') || strcmp(pm,'none')
                obj.pixelmode = pm;
            end
        end
    end
    
    methods (Access=public)
        % N-D phase correlation
        % F-> moving, G-> fixed
        function [r,Q] = FourierShiftND(obj,F,G,dims,pixmode,averaging)
            if nargin < 6
                averaging = obj.averaging;
            end
            
            if nargin < 5
                pixmode   = obj.pixelmode;
            end
            
            if nargin < 4
                dims = 1:ndims(F);
            end
            
            alpha = 0.15; % Tukey alpha
            
            ndim  = numel(dims);
            sz    = size(F);
            
            f = single(F);
            g = single(G);
            
            % Calculate a Tukey window for tapering edges
            for i = 1:ndim
                f = f.*single(obj.TukeyWindow(sz(dims(i)),alpha,dims(i)));
                g = g.*single(obj.TukeyWindow(sz(dims(i)),alpha,dims(i)));
            end
            
            for i = 1:ndim
                f = fft(f,[],dims(i));
                g = fft(g,[],dims(i));
            end
            
            % Calculate cross power spectrum (Reddy et al, IEEE 1996)
            q = (conj(f).*g)./max(eps,abs(f.*g));
            
            % Inverse fft of cross power spectrum
            Q = q;
            for i = 1:ndim
                Q = ifft(Q,[],dims(i));
            end
            
            Q = abs(Q);
            
            for i = 1:ndim
                Q = fftshift(Q,dims(i)); % Shift zero freq to center before fitting
            end
            
            if averaging
                sdims       = 1:3;
                sdims(dims) = [];
                if ~isempty(sdims)
                    Q       = sum(Q,sdims);
                end
            end
            
            [~,midx]   = max(Q,[],'all','linear');
            [ri,ci,zi] = ind2sub(size(Q),midx);
            wdim       = [ri,ci,zi];
            
            switch lower(pixmode)
                case 'subpixel'
                    wdw             = 7;    % Fitting window in pixels (+/-)
                    idx             = {[wdim(1)],[wdim(2)],[wdim(3)]};
                    for i = dims
                        rl          = BoundCheck([wdim(dims(i))-wdw wdim(dims(i))+wdw],sz(dims(i)));
                        idx{dims(i)} = rl(1):rl(2);
                        cpix(i)     = rl(1);
                    end
                    
                    grab            = Q(idx{:});
                    
                    switch ndim
                        case 1
                            [fparams,sparams] = obj.FitGauss1D(grab,cpix,'lsq');
                            if ~isempty(fparams)
                                r = fparams(1,2)-sz(dims)/2-1+mod(sz(dims),2)*0.5;
                            else
                                r = sparams(1,2)-sz(dims)/2-1+mod(sz(dims),2)*0.5;
                                obj.AddMessage('Fit unsuccessful: Reverted to centroid.',1);
                            end
                        case 2
                            [fparams,sparams] = obj.FitGauss2D(grab,cpix,'lsq');
                            if ~isempty(fparams)
                                r(1) = fparams(1,3)-sz(dims(1))/2-1+mod(sz(dims(1)),2)*0.5;
                                r(2) = fparams(1,2)-sz(dims(2))/2-1+mod(sz(dims(2)),2)*0.5;
                            else
                                r(1) = sparams(1,3)-sz(dims(1))/2-1+mod(sz(dims(1)),2)*0.5;
                                r(2) = sparams(1,2)-sz(dims(2))/2-1+mod(sz(dims(2)),2)*0.5;
                                obj.AddMessage('Fit unsuccessful: Reverted to centroid.',1);
                            end
                        case 3
                            obj.AddMessage('3D fitting not implemented: Reverted to pixel precision.',1)
                            r = wdim(dims)-sz(dims)/2-1+mod(sz(dims),2)*0.5;
                    end
                case 'pixel'
                    r = wdim(dims)-sz(dims)/2-1+mod(sz(dims),2)*0.5;
                case 'none'
                    r = [];
            end
            
            % Bound checks for grab box
            function rl = BoundCheck(rl,sz)
                if rl(1)<1;  rl = rl+abs(min(rl))+1; end
                if rl(2)>sz; rl = rl+(sz-rl(2));     end                
            end
        end
        
        % Fit to 1D gaussian model
        function [fparams,sparams,imfit] = FitGauss1D(obj,grabs,cpix,fitmethod)
            if nargin<4 % Set default method
                fitmethod = 'lsq'; 
            end
            
            if nargin<3 % Set default corner pixels to null
                cpix = [];
            end
            
            ti=tic;
            warning off
            
            fparams = zeros(size(grabs,3),4); % Preallocate parameter matrices
            
            % Set up coordinate grid for fitting
            X = 1:length(grabs);
            
            if ~strcmpi(fitmethod,'nlinfit') && ~strcmpi(fitmethod,'lsq')
                error('Unknown method. Must be nlinfit or lsq.');
            end
            
            % Set fit options
            switch lower(fitmethod)
                case 'nlinfit'
                    opts = statset('nlinfit');
                case 'lsq'
                    opts = optimoptions('lsqcurvefit');
                    opts.Algorithm     = 'levenberg-marquardt';
                    opts.ScaleProblem  = 'jacobian';
                    opts.FunctionTolerance = 1e-6;
                    opts.StepTolerance = 1e-6;
                    opts.Display       = 'none';
                    lbound = [0,0,0,0];
                    ubound = [inf,size(grabs,2),inf,inf];
            end
            
            for i=1:size(grabs,3)
                testgrab = double(grabs(:,:,i));
                
                % Estimate starting parameters
                xc = sum(X.*testgrab,'all')./sum(testgrab,'all');
                xm = max(testgrab,[],'all');
                xw = find((xm>0.5*max(xm))==1);
                xw = max(xw)-min(xw)+1; % Estimate x FWHM
                xw = xw/sqrt(8*log(2)); % Convert from FWHM to sigma
                
                sparams(i,:) = [max(testgrab,[],'all'),xc,xw,min(testgrab(:))];
                
                try
                    switch lower(fitmethod)
                        case 'nlinfit'
                            lastwarn('','');
                            fparams(i,:) = nlinfit(X,testgrab(:),@gaussian1Nlin,sparams,opts);
                            [~,warnId]   = lastwarn();
                            if ~isempty(warnId) % If the fit didn't converge, flag it for removal
                                fparams(i,:) = NaN;
                            end
                        case 'lsq'
                            [fparams(i,:),~,~,exitflag] = lsqcurvefit(@gaussian1,sparams,X,testgrab,lbound,ubound,opts);
                            if exitflag < 1 || sum(fparams(i,:)==(ubound+eps)) || sum(fparams(i,:)==(lbound-eps)) % If the fit didn't converge or is outside the boundaries, flag it for removal
                                fparams(i,:) = NaN;
                            end
                    end
                    
                    if nargout>2
                        imfit(i,:) = gaussian1(fparams(i,:),X);
                    end
                catch ME
                    fparams(i,:) = NaN;
                    warning on
                    warning(ME.message)
                end
                
                if ~isempty(cpix) % Adjust to coordinate frame of original image if those coordinates are available
                    fparams(i,2) = fparams(i,2)+cpix(i,1)-1;
                    sparams(i,2) = sparams(i,2)+cpix(i,1)-1;                        
                end
            end
            
            fparams = fparams(sum(isnan(fparams),2)==0,:); % Clear out non-convergent fits
            
            obj.AddMessage([num2str(size(fparams,1)),' psfs successfully fitted in ',num2str(toc(ti)),'s.'])
            
            warning on
            
            %% Gaussian fitting functions
            % Lsqcurvefit (optimization toolbox)
            % ---------------------------------
            % Gaussian 1D (1d pdf)
            function f = gaussian1(b,xy)
                f = b(1)*exp(-(((xy(:,:,1)-b(2)).^2)./(2*b(3).^2)))+b(4);
            end
            
            % Gaussian 1D (Gaussian beam equation)
            function f = gaussian1w(b,xy)
                f = b(1)*exp(-2*(((xy(:,:,1)-b(2)).^2)./(b(3).^2)))+b(4);
            end
            
            % Nlinfit (stats toolbox)
            % ---------------------------------
            % Gaussian 1D (1d pdf)
            function f = gaussian1Nlin(b,xy)
                f = b(1)*exp(-(((xy(:,:,1)-b(2)).^2)./(2*b(3).^2)))+b(4);
                f = f(:);
            end
            
            % Gaussian 1D (Gaussian beam equation)
            function f = gaussian1wNlin(b,xy)
                f = b(1)*exp(-2*(((xy(:,:,1)-b(2)).^2)./(b(3).^2)))+b(4);
                f = f(:);
            end
        end
        
        % Fit to 2D gaussian model
        % -Maps to image coordinate system if the top left corner of the grab windows (cpix) are passed into the function
        function [fparams,sparams,imfit] = FitGauss2D(obj,grabs,cpix,fitmethod)
            if nargin<4 % Set default method
                fitmethod = 'lsq'; 
            end
            
            if nargin<3 % Set default corner pixels to null
                cpix = [];
            end
            
            ti=tic;
            warning off
            
            fparams = zeros(size(grabs,3),7); % Preallocate parameter matrices
            sparams = zeros(size(grabs,3),7);
            
            % Set up coordinate grid for fitting
            [X,Y]     = meshgrid(1:size(grabs,2),1:size(grabs,1));
            XY(:,:,1) = X;
            XY(:,:,2) = Y;
            
            if ~strcmpi(fitmethod,'nlinfit') && ~strcmpi(fitmethod,'lsq')
                error('Unknown method. Must be nlinfit or lsq.');
            end
            
            % Set fit options
            switch lower(fitmethod)
                case 'nlinfit'
                    opts = statset('nlinfit');
                case 'lsq'
                    opts = optimoptions('lsqcurvefit');
                    opts.Algorithm     = 'levenberg-marquardt';
                    opts.ScaleProblem  = 'jacobian';
                    opts.FunctionTolerance = 1e-6;
                    opts.StepTolerance = 1e-6;
                    opts.Display       = 'none';
                    lbound = [0,0,0,0,0,0,-inf];
                    ubound = [inf,size(grabs,2),size(grabs,1),inf,inf,inf,inf];
            end
            
            for i=1:size(grabs,3)
                testgrab = double(grabs(:,:,i));
                
                % Estimate starting parameters
                yc = sum(Y.*testgrab,'all')./sum(testgrab,'all');
                xc = sum(X.*testgrab,'all')./sum(testgrab,'all');
                % [yc,xc] = find(testgrab==max(testgrab,[],'all'),1); % Estimate start coordinate
                ym = max(testgrab,[],1);
                xm = max(testgrab,[],2);
                yw = find((ym>0.5*max(ym))==1);
                yw = max(yw)-min(yw)+1; % Estimate y FWHM
                yw = yw/sqrt(8*log(2)); % Convert from FWHM to sigma
                xw = find((xm>0.5*max(xm))==1);
                xw = max(xw)-min(xw)+1; % Estimate x FWHM
                xw = xw/sqrt(8*log(2)); % Convert from FWHM to sigma
                
                sparams(i,:) = [max(testgrab,[],'all'),xc,yc,xw,yw,min(testgrab(:)),0];
                
                try
                    switch lower(fitmethod)
                        case 'nlinfit'
                            lastwarn('','');
                            fparams(i,:) = nlinfit(XY,testgrab(:),@gaussian2bvNlin,sparams,opts);
                            [~,warnId]   = lastwarn();
                            if ~isempty(warnId) % If the fit didn't converge, flag it for removal
                                fparams(i,:) = NaN;
                            end
                        case 'lsq'
                            [fparams(i,:),~,~,exitflag] = lsqcurvefit(@gaussian2bv,sparams,XY,testgrab,lbound,ubound,opts);
                            if exitflag < 1 || sum(fparams(i,:)==(ubound+eps)) || sum(fparams(i,:)==(lbound-eps)) % If the fit didn't converge or is outside the boundaries, flag it for removal
                                fparams(i,:) = NaN;
                            end
                    end
                    
                    if nargout>2 
                        imfit(:,:,i) = gaussian2bv(fparams(i,:),XY);
                    end
                catch ME
                    fparams(i,:) = NaN;
                    warning on
                    warning(ME.message)
                end
                
                if ~isempty(cpix) % Adjust to coordinate frame of original image if those coordinates are available
                    fparams(i,2) = fparams(i,2)+cpix(i,2)-1;
                    fparams(i,3) = fparams(i,3)+cpix(i,1)-1;
                    sparams(i,2) = sparams(i,2)+cpix(i,2)-1;
                    sparams(i,3) = sparams(i,3)+cpix(i,1)-1;                        
                end
            end
            
            fparams = fparams(sum(isnan(fparams),2)==0,:); % Clear out non-convergent fits
            
            obj.AddMessage([num2str(size(fparams,1)),' psfs successfully fitted in ',num2str(toc(ti)),'s.'])
            
            warning on
            
            %% Gaussian fitting functions
            % Lsqcurvefit (optimization toolbox)
            % ---------------------------------
            % Gaussian 2D (bivariate)
            function f = gaussian2bv(b,xy)
                x   = (xy(:,:,1)-b(2)).^2./b(4).^2;
                y   = (xy(:,:,2)-b(3)).^2./b(5).^2;
                xyt = 2*b(7).*(xy(:,:,1)-b(2)).*(xy(:,:,2)-b(3))./(b(4).*b(5));
                
                z   = x+y-xyt;
                
                f   = b(1)*exp((-1/2*z)./(1-b(7).^2))+b(6);
            end
            
            % Gaussian 2D (2d pdf)
            function f = gaussian2(b,xy)
                f = b(1)*exp(-(((xy(:,:,1)-b(2)).^2)./(2*b(4).^2)+((xy(:,:,2)-b(3)).^2)./(2*b(5).^2)))+b(6);
            end
            
            % Gaussian 2D (Gaussian beam equation)
            function f = gaussian2w(b,xy)
                f = b(1)*exp(-2*(((xy(:,:,1)-b(2)).^2)./(b(4).^2)+((xy(:,:,2)-b(3)).^2)./(b(5).^2)))+b(6);
            end
            
            % Nlinfit (stats toolbox)
            % ---------------------------------
            % Gaussian 2D (bivariate)
            function f = gaussian2bvNlin(b,xy)
                x   = (xy(:,:,1)-b(2)).^2./b(4).^2;
                y   = (xy(:,:,2)-b(3)).^2./b(5).^2;
                xyt = 2*b(7).*(xy(:,:,1)-b(2)).*(xy(:,:,2)-b(3))./(b(4).*b(5));
                
                z   = x+y-xyt;
                
                f   = b(1)*exp((-1/2*z)./(1-b(7).^2))+b(6);
                f   = f(:);
            end
            
            % Gaussian 2D (2d pdf)
            function f = gaussian2Nlin(b,xy)
                f = b(1)*exp(-(((xy(:,:,1)-b(2)).^2)./(2*b(4).^2)+((xy(:,:,2)-b(3)).^2)./(2*b(5).^2)))+b(6);
                f = f(:);
            end
            
            % Gaussian 2D (Gaussian beam equation)
            function f = gaussian2wNlin(b,xy)
                f = b(1)*exp(-2*(((xy(:,:,1)-b(2)).^2)./(b(4).^2)+((xy(:,:,2)-b(3)).^2)./(b(5).^2)))+b(6);
                f = f(:);
            end
        end
        
        % Generate a Tukey window of sample N and parameter alpha (for fft pre-processing)
        function w = TukeyWindow(~,N,alpha,dim)
            k          = round(alpha*N/2);
            w          = ones(N,1);
            w(1:k)     = 0.5*(1-cos(2*pi*linspace(0,k,k)/(alpha*N+1)));
            w(N-k+1:N) = w(k:-1:1);
            
            if dim > 1
                rsdim      = zeros(1,dim)+1;
                rsdim(dim) = N;
                w          = reshape(w,rsdim);
            end
        end
    end
    
    methods (Access=private)
        % Add message to command window (if warn is 1 or 2, add message and ignore verbose)
        function AddMessage(obj,msg,warn)
            if nargin < 3
                warn = 0;
            end
            
            switch warn
                case 0
                    if obj.verbose
                        disp([class(obj),': ',msg])
                    end
                case 1 % Always display warnings
                    warning on
                    warning([class(obj),': Warning: ',msg])
                case 2
                    error([class(obj),': Error: ',msg])
            end
        end
    end
end