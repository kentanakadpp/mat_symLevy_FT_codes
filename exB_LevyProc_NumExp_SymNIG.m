%% Function mu
mu  = @(x) x .* besselk(1,x)/pi;

%% The (half) Fourier transform of mu
% FTmu = @(x) (pi - 2*1i*(x*sqrt(1+x^2)+asinh(x)))/(2*pi*(1+x^2)^(3/2));
%% The "single" indefinite integral of FTmu
% SIFTmu = @(x) (0.5 - 1i*asinh(x)/pi)*x/sqrt(1+x^2);
%% The "double" indefinite integral of FTmu
% DIFTmu = @(x) (sqrt(1+x^2)-1)/2 + 1i*(x - sqrt(1+x^2)*asinh(x))/pi;

%% The exact solution 
T = 3; % Time
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
i_frst = 7;
i_last = 12;
M_list = zeros(1,i_last-i_frst+1);
ErrArr = zeros(1,i_last-i_frst+1);
ErrPrt = zeros(1,i_last-i_frst+1);
CpTime = zeros(1,i_last-i_frst+1);
for i=i_frst:i_last
    %% Parameters depending on M
    % (for DE-FT)
    M  = 2^i;
    Mm = M/2;
    h  = log(10^3*M)/M; % 0.025;
    % (for the first indefinite integral)
    L = floor(M/4);
    % (for Continuous Euler FT)
    N = floor(L/2);
        
    %% Numerical solution of the symmetric Levy PIDE
    tic;
    % (including two-times indefinite integration)
    NumSol = LevyPDE_SolFunc_IItwo(mu, T, M, Mm, h, beta, eps, Mmax, wd, wu, d, L, N); 
    CpTime(i-i_frst+1) = toc;
    
    %% Outputs
    % Exact solution
    h2 = wu/N;
    ES = zeros(1,2*N);
    for m=1:2*N
        ES(m) = ESol((m-N)*h2);
    end
    % M
    M_list(i-i_frst+1) = M;
    % Max Log Error
    Log10Err = log10(abs(NumSol - ES));
    ErrArr(i-i_frst+1) = max(Log10Err);
    ErrPrt(i-i_frst+1) = max(horzcat(Log10Err(1:1+floor((wu-wd)/h2)),Log10Err(2*N-floor((wu-wd)/h2):2*N)));
end

%% Error plot
hold on;
plot(M_list, ErrArr, '--sb', 'MarkerSize',8, 'MarkerFaceColor','b', 'LineWidth',3);
plot(M_list, ErrPrt, '--<r', 'MarkerSize',8, 'MarkerFaceColor','r', 'LineWidth',3);
ylim([-15,-1]);
legend(['error on [-',num2str(wu),',',num2str(wu),']'],...
    ['error on [-',num2str(wu),',-',num2str(wd),'] \cup [',num2str(wd),',',num2str(wu),']'],...
    'Location','SouthWest');
grid on;
set(gca,'FontName','Times','FontSize',20,'FontWeight','bold'); 
xlabel('M');
ylabel('log_{10}(Max Error)');
title(['Errors for the NIG process (t = ',num2str(T),')']);
commstr = strcat('print -depsc fig_NIG_error_t', num2str(T), '.eps');
eval(commstr);
hold off;

%% Time plot
plot(M_list, CpTime, '--og', 'MarkerSize',8, 'MarkerFaceColor','g', 'LineWidth',3);
ylim([0,1]);
grid on;
set(gca,'FontName','Times','FontSize',20,'FontWeight','bold'); 
xlabel('M');
ylabel('Time (sec)');
title(['Computational times for the NIG process (t = ',num2str(T),')']);
commstr = strcat('print -depsc fig_NIG_cptime_t', num2str(T), '.eps');
eval(commstr);
hold off;
