function Q = FourierShift2D(F,G)
        F = single(F);
        G = single(G);
        
        sz = size(F);
        
        % Calculate a Tukey window for tapering edges
        w = single(TukeyWindow(sz(1),0.1)*TukeyWindow(sz(2),0.1)');
        
        % Apply Tukey window
        f = F.*w;
        g = G.*w;
        
        f = fft2(f);
        g = fft2(g);
        
        % Calculate cross power spectrum (Reddy et al, IEEE 1996)
        q = (conj(f).*g)./max(eps,abs(f.*g));
        
        % Inverse fft of cross power spectrum
        Q = ifft2(q);
        
        Q = abs(Q);
        Q = fftshift(Q); % Shift zero freq to center before fitting
        
        % Generate a Tukey window of sample N and parameter alpha (for fft pre-processing)
        function w = TukeyWindow(N,alpha)
            k          = round(alpha*N/2);
            w          = ones(N,1);
            w(1:k)     = 0.5*(1-cos(2*pi*linspace(0,k,k)/(alpha*N+1)));
            w(N-k+1:N) = w(k:-1:1);
        end
    end