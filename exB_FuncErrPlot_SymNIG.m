%% Function mu
mu  = @(x) x .* besselk(1,x)/pi;

%% The (half) Fourier transform of mu
% FTmu = @(x) (pi - 2*1i*(x*sqrt(1+x^2)+asinh(x)))/(2*pi*(1+x^2)^(3/2));
%% The "single" indefinite integral of FTmu
% SIFTmu = @(x) (0.5 - 1i*asinh(x)/pi)*x/sqrt(1+x^2);
%% The "double" indefinite integral of FTmu
% DIFTmu = @(x) (sqrt(1+x^2)-1)/2 + 1i*(x - sqrt(1+x^2)*asinh(x))/pi;

for t=1:3
    %% The exact solution 
    T = t; % Time
    ESol = @(x)  T * exp(T) * besselk(1,sqrt(T^2+x^2))/(pi*sqrt(T^2+x^2));

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
    % (for the first indefinite integral)
    L = floor(M/4);
    % (for Continuous Euler FT)
    N = floor(L/2);
        
    %% Numerical solution of the symmetric Levy PIDE
    % (including two-times indefinite integration)
    NumSol = LevyPDE_SolFunc_IItwo(mu, T, M, Mm, h, beta, eps, Mmax, wd, wu, d, L, N); 
    
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
title(['Solution for the NIG process']);
commstr = strcat('print -depsc fig_NIG_exact.eps');
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
ylabel('log_{10}(Max Error)');
title(['Errors for the NIG process (M = ', num2str(M),')']);
commstr = strcat('print -depsc fig_NIG_error_xeplane.eps');
eval(commstr);
hold off;
