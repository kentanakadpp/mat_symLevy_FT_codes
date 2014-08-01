%% Function mu
mu  = @(x) exp(-x);

%% The (half) Fourier transform of mu
% FTmu = @(x) 1/(1+1i*x);
%% The indefinite integral of FTmu
% IFTmu = @(x) -1i*log(1+1i*x);

for t=1:3
    %% The exact solution 
    T = t; % Time
    ESol = @(x) 2^(0.5-T)*abs(x)^(T-0.5)*besselk(0.5-T,abs(x))/(sqrt(pi)*gamma(T));

    %% ========================================================================

    %% Parameters
    % (for DE-FT)
    beta = 0.25;
    % (for Nonuniform FT)
    eps = 1e-10;
    Mmax = 20;
    % (for Continuous Euler FT)
    wd = 2;
    wu = 5;
    d = 1;

    %% Numerical solutions of the Levy PIDE and their errors

    %% Parameters depending on M

    % (for DE-FT)
    M  = 2^11;

    Mm = M/2;
    h  = log(10^3*M)/M; % 0.025;
    % (for Continuous Euler FT)
    N = floor(M/4);
        
    %% Numerical solution of the symmetric Levy PIDE
    NumSol = LevyPDE_SolFunc_IIone(mu, T, M, Mm, h, beta, eps, Mmax, wd, wu, d, N);
    
    %% Outputs
    % Exact solution
    h2 = wu/N;
    x = [-N+1:N]*h2;
    ES = zeros(1,2*N);
    for m=1:2*N
        ES(m) = ESol((m-N)*h2);
    end
    % Max Log Error
    Log10Err = log10(abs(NumSol - ES));
    
    if(t==1)
       ES1 = ES;
       LErr1 = Log10Err;
    elseif(t==2)
       ES2 = ES;
       LErr2 = Log10Err;
    elseif(t==3)
       ES3 = ES;
       LErr3 = Log10Err;
    end    
    
end

%% Exact solution plot
hold on;
plot(x, ES1, '-k', 'LineWidth',3);
plot(x, ES2, '--b', 'LineWidth',3);
plot(x, ES3, '-.r', 'LineWidth',3);
ylim([0,0.6]);
legend('t = 1', 't = 2', 't = 3');
grid on;
set(gca,'FontName','Times','FontSize',20,'FontWeight','bold'); 
xlabel('x');
ylabel('p(x,t)');
title(['Solution for the VG process']);
commstr = strcat('print -depsc fig_VG_exact.eps');
eval(commstr);
hold off;
    
%% Error plot
hold on;
plot(x, LErr1, '-k', 'LineWidth',3);
plot(x, LErr2, '--b', 'LineWidth',3);
plot(x, LErr3, '-.r', 'LineWidth',3);
ylim([-15,-1]);
legend('t = 1', 't = 2', 't = 3');
grid on;
set(gca,'FontName','Times','FontSize',20,'FontWeight','bold'); 
xlabel('x');
ylabel('log_{10}(Error)');
title(['Errors for the VG process (M = ', num2str(M),')']);
commstr = strcat('print -depsc fig_VG_error_xeplane.eps');
eval(commstr);
hold off;
