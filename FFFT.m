%%%%% [Fractional FFT] %%%%%
% Routine to compute the sums:
% \sum_{n=1}^{2N} \alpha(n) \exp(i 2\pi m (n-N) a)
% (m = -N+1,...,N).

function FT = FFFT(alpha, a, N)

e1 = zeros(1, 2*N);
e2 = zeros(1, 2*N);
y1 = zeros(1, 2*N);
y2 = zeros(1, 2*N);
z1 = zeros(1, 2*N);
z2 = zeros(1, 2*N);
for n=1:2*N
    e1(n) = exp(-1i*2*pi*(n-N)*(N-1)*a + 1i*pi*(n-1)^2*a);
    e2(n) = exp(-1i*2*pi*(n-1)*(N-1)*a + 1i*pi*(n-1)^2*a);
    y1(n) = alpha(n) * e1(n);
    z1(n) = exp(-1i*pi*(n-1)^2*a);
    z2(n) = exp(-1i*pi*(2*N-(n-1))^2*a);
end
y = horzcat(y1,y2);
z = horzcat(z1,z2);
Y = fft(y,4*N);
Z = fft(z,4*N);
X = Y .* Z;
x = ifft(X,4*N);

FT = e2 .* x(1:2*N);

end
