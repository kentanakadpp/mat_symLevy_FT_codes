%%%%% [DE-FT]+[NFFT] %%%%%
% Routine to compute the (half) Fourier transform 
% using the double-exponential formula by Ooura 
% combined with non-uniform FFT.
% 
% [Input]
% f = function to be transformed
% [Output function values]
% Ff = the (half) Fourier transform of f for (n-1)*ht, where n = 1,...,M. 

function Ff = DE_NFFT(f, z0, M, Mm, h, beta, ht, eps, Mmax)

%% Making DE coefficients
alpha = beta/(1+log(1+pi/(z0*h))/(4*z0*h));
x = zeros(1,M);
PhiMu = zeros(1,M);
for j=1:M
    [x(j),dj,yj]=phi((j-Mm)*h, alpha, beta);
    PhiMu(j)=-2*pi*1i/z0 * f(pi/(z0*h)*x(j)) * dj * sin(pi/(2*h)*yj) * exp(pi*1i/(2*h)*yj);
end
t = (M*ht/(2*z0*h))*x;

% tic;
%% Making index set
j_max = ones(1,M);
j_min = ones(1,M);
Mlmt = M;
for m=1:M
    if(m-1 > ceil(t(M))+2*Mmax)
        Mlmt = m;
        break;
    end
    % j_max(m)
    if(m==1)
        j_max_init = 1;
    else
        j_max_init = j_max(m-1);
    end
    for j=j_max_init:M
        if(m-1 >= floor(t(j)))
            j_max(m) = j;
        else
            break;
        end
    end
    % j_min(m)
    if(m==1)
        j_min_init = 1;
    else
        j_min_init = j_min(m-1);
    end
    for j=j_min_init:M
        if(m-1 >= ceil(t(j))+2*Mmax)
            j_min(m) = j;
        else
            break;
        end
    end
end
% toc;

%% Gaussian interpolation
tau = -log(eps)/(pi^2);
hc = 1; % step width for discritization of the Fourier integral of the Gaussian kernel
Msft = floor(M/4); % shift of the indices of GdgMu() to decrease Ead2()
GgdMu = zeros(1,M);
for m = 1:Mlmt
    mid_j = floor((j_min(m)+j_max(m))/2);
    tmp_Ggd_l = 0;
    for j = j_min(m):mid_j
        tmp_Ggd_l = tmp_Ggd_l + (PhiMu(j)*hc/(2*pi)) * exp(-1i*2*pi*Msft*t(j)/M) * exp(- hc^2 * (m-1-Mmax - t(j))^2/(4*tau));
    end
    tmp_Ggd_u = 0;
    for j = j_max(m):-1:(mid_j+1)
        tmp_Ggd_u = tmp_Ggd_u + (PhiMu(j)*hc/(2*pi)) * exp(-1i*2*pi*Msft*t(j)/M) * exp(- hc^2 * (m-1-Mmax - t(j))^2/(4*tau));
    end    
    GgdMu(m) = (tmp_Ggd_l + tmp_Ggd_u) * exp(1i*2*pi*Msft*(m-1)/M);
end

%% Adjustment coefficients
Ead1 = zeros(1,M);
Ead2 = zeros(1,M);
for m=1:M
    Ead1(m) = exp(1i*2*pi*Mmax*(m-1-Msft)/M);
    Ead2(m) = sqrt(pi/tau)*exp(tau*(2*pi*(m-1-Msft)/(M*hc))^2);
end

%% FFT
Ff = Ead2 .* Ead1 .* fft(GgdMu);

end
