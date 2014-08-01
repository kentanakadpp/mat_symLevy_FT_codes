%%%%% [SG indefinite integral] %%%%%
% Routine to compute the indefinite integral of an even function
%  (a symmetric function w.r.t. the origin) 
% using Sinc-Gauss sampling formula.
% 
% [Input function values]
% IntegrandPlus = array of the integrand values at (n-1)*ht, where n = 1,...,2N+1.
% [Output function values]
% SGII = array of the integral values at (n-1)*ht, where n = 2,...,N+1.

function SGII = SG_IndefInt_sym(IntegrandPlus, N, ht, Nn) % Nn > 2*N

%% Generating function of G's
r = sqrt(N/pi);    
Fsg = @(w) 0.5 * ( erf(r*(w+pi)/sqrt(2)) - erf(r*(w-pi)/sqrt(2)) );

%% Computation of G's using fractional FFT
ha = 2*pi/Nn; % ha = 2*pi/Nn
a = ha/(2*pi);
F = zeros(1, 2*Nn);
for m=1:2*Nn
    n = m-Nn;
    if(n==0)
        tmp_sinc = 1.0;
    else
        tmp_sinc = sin(n*ha/2)/(n*ha/2);
    end
    F(m) = (ha/(2*pi)) * Fsg(n*ha) * tmp_sinc * exp(-1i*n*ha/2);
end

g = FFFT(F, a, Nn); % Fractional FFT

G = zeros(1,2*N+1);
for i=1:N
   G(N+i+1) = G(N+i) + g(Nn+i);
   G(N-i+1) = - G(N+i+1);
end

%% Symmetrization of the integrand array
SymF = zeros(1,4*N);
for l=1:2*N-1
    SymF(l) = conj(IntegrandPlus(2*N+1-l));
end
for l=2*N:4*N
    SymF(l) = IntegrandPlus(l+1-2*N);
end

%% Sinc-Gauss indefinte integration for an even function
H = zeros(1,N);
H(1)=0.0;
for l=2:N
    H(l) = H(l-1) + ht * (-G(1) * SymF(3*N+l-1) + G(2*N+1) * SymF(N+l-1));
end

IP = 0.0;
for k=-N+1:N
    IP = IP - ht * SymF(2*N+k) * G(N+1-k);    
end
I = ones(1,N) * IP;

% Convolution by FFT
esft = zeros(1,4*N);
mu = zeros(1,4*N);
gr = zeros(1,4*N);
for k=-2*N+1:2*N
    esft(2*N+k) = exp(1i*pi*(2*N-1)*(2*N+k-1)/(2*N));
    mu(2*N+k) = SymF(2*N+k) * esft(2*N+k);
    if((k<=-N)||(N+1<=k))
        gr(2*N+k) = 0.0;
    else
        gr(2*N+k) = G(N+1+k) * esft(2*N+k);
    end
end
ecnt = exp(-1i*pi*(2*N-1)^2/(2*N));
DFmu = ecnt * (esft .* fft(mu,4*N)); % FFT
DFgr = ecnt * (esft .* fft(gr,4*N)); % FFT
prod = DFmu .* DFgr .* conj(esft);
IDFmg = conj(ecnt) * (conj(esft) .* ifft(prod,4*N)); % IFFT
IIF = ht * IDFmg(1,2*N+1:3*N);

SGII = IIF + I + H;

end
