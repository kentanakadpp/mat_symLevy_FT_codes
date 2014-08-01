%%%%% [Continuous Euler + Fractinal FFT] %%%%%
% Routine to compute the sums:
% \sum_{n=1}^{2N} (h1/(2*pi)) * w((n-N)*h1) * f(n) \exp(i m (n-N) h1 h2)
% (m = -N+1,...,N),
% where
% h1 = ht, 
% h2 = wu/N,
% f = array of function values at (n-N)*ht, where n = 1,...,2N

function CEFT = ContEulerFFFT(f, wd, wu, N, ht)

%% Parameters and weight function
h1 = ht;
h2 = wu/N;
p = sqrt(N*h1/wd);
q = sqrt(wd*N*h1/4.0);
w = @(x) 0.5*erfc(abs(x)/p-q);

%% Fractional FFT
a = h1*h2/(2*pi);
F = zeros(1, 2*N);
for n=1:2*N
    F(n) = (h1/(2*pi)) * w((n-N)*h1) * f(n);
end

CEFT = FFFT(F, a, N); % Fractional FFT

end
