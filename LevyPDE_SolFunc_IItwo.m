function NumSol= LevyPDE_SolFunc_IItwo(mu, T, M, Mm, h, beta, eps, Mmax, wd, wu, d, L, N)
    %% STEP 1: [DE-FT]+[NFFT]         
    
    % Parameters for Continuous Euler FT (Needed here to determine "ht")
    ht = sqrt(2*pi*d*(wd+wu)/(wd^2*N));

    % DE-FT + NFFT
    z0 = M*ht/30;
    z1 = M*ht/3.6;
    NumFTmu0 = DE_NFFT(mu, z0, M, Mm, h, beta, ht, eps, Mmax);
    NumFTmu1 = DE_NFFT(mu, z1, M, Mm, h, beta, ht, eps, Mmax);
    NumFTmu = horzcat(NumFTmu0(1:floor(M/16)),NumFTmu1(floor(M/16)+1:M));

    %% STEP 2: [SG indefinite integral]
    Nn = M; % Nn > 2*N
    IIF = SG_IndefInt_sym(NumFTmu, L, ht, Nn); % different from 'IIone'
    IIF = horzcat([0], IIF);                   % different from 'IIone'
    IIF = -1i * SG_IndefInt_sym(1i * IIF, N, ht, Nn);     % different from 'IIone'

    %% STEP 3: [Continuous Euler FT]
    % Exponential of the indefinite integral multiplied by T and its symmetrization
    ExpChr = exp(-T*2*real(IIF));              % different from 'IIone'

    SymEC = zeros(1,2*N);
    for n=1:N-1
        SymEC(n) = ExpChr(N-n);
    end
    SymEC(N) = 1.0;
    for n=N+1:2*N
        SymEC(n) = ExpChr(n-N);
    end

    % Continuous Euler + Fractional FFT
    NumSol = ContEulerFFFT(SymEC, wd, wu, N, ht);

end
