close all; clear all; clc;

%==========================================================================
% Solves the income fluctuation model 
%
% Last updated 11/7/2013
%==========================================================================

global rho pmin draw_eta draw_eps draw_u B beta mu R nwgrid xigrid cmax;

% Initialize the random number generator
s = RandStream.create('mt19937ar','seed',5489);

%Stopping criterion. Stop when 
% max_node|u'(c^N(node)-c^{N-1}(node)| < Tol
Tol = 1e-4;

%Parameters
beta   = 0.95;                  % discount rate
mu     = 2;                     % CRRA risk aversion
r      = 0.02;                  % interest rate
R      = 1 + r;                 

exercise = 1;                   % choose parameter values

if exercise == 1;               % baseline exercise
    mean_eta = -5;          % expected minimum income close to zero
    B = 0;                      % no borrowing is close to the natural debt limit
    rho_vec = [1 0.95 0.5];     % persistence of income shock
elseif exercise == 2;           % tightening the borrowing limit relative to the baseline
    mean_eta = -5;          
    B = -10;                    % precommitted spending
    rho_vec = [1 0.95];
else
    mean_eta = 0;           % raising the expected minimum level of income relative to the baseline
    B = 0;                      
    rho_vec = [1 0.95];
end;    

% Simulate points for Monte Carlo integration
% Income process: 
% Y_t = xi_t * eps_t + eta_t
% eta_t ~ LN(mean_eta,sig_eta^2)
% ln xi_t ~ rho*ln(xi_{t-1}) + ln u_t, u_t ~LN (0,sig_xi^2)
% eps_t ~ LN(0,sig_eps^2)
N = 1000;                        % number of draws for Monte Carlo Simulation 
sig_eps = 0.05;                 % std of multiplicative transitory shock
sig_u = 0.01;                  % std of persistent component
sig_xi0 = 0.15;                 % xi0 ~LN(0, sig_xi0^2)
sig_eta = 0.001;                % std of additive transitory shockmean_eta = 0; 
draws = randn(N,3);             % random draws for each shock
draw_eta = exp(mean_eta + draws(:,1)*sig_eta); % simulated additive transitory shock  
draw_eps = exp(draws(:,2)*sig_eps); % simulated multiplicative transitory shock  
draw_u = exp(draws(:,3)*sig_u);   % simulated innovation to the persistent component

for rr=1:length(rho_vec);
    rho = rho_vec(rr);
    
    %Create grid for xi 
    ximax   = 25;                  %maximum of grid  
    ximin   = 1e-4;                 %minium 
    xisize  = 50;                   %grid size
    xiscale = 2;                    %grid scale, increasing the scale concentrate points in lower region of the grid
    xigrid =  makegrid(ximin, ximax, xisize, xiscale);   %makegrid is a function that creates the grid

    %Create networth grid
    nwmax   = 100;                  %maximum of grid  
    nwmin   = -B*R + ximin;         %minium is the borrowing limit
    nwsize  = 100;                   %grid size
    nwscale = 2;                    %grid scale, increasing the scale concentrate points in lower region of the grid
    nwgrid =  makegrid(nwmin, nwmax, nwsize, nwscale);   %makegrid is a function that creates the grid

    %Construct an initial guess as a fraction of the maximum feasible consumption 
    cmax = nwgrid + B; %maximum feasible consumption
    cmax = repmat(cmax,1,xisize);
    cfrac = 0.1;                        
    c0 = cfrac * cmax; % consume 10% of income
    p0 = u1(c0,mu); % price function 
    pmin  = u1(cmax,mu); % marginal utility at the maximum consumption level
    
    %name to save result
    fname = ['Exercise_', num2str(exercise), '_rho',num2str(round(100*rho))];
   
    %policy iteration to find the optimal consumption policy
    tic;    
    err = 1;count = 1; maxiter = 1000;

    while err > Tol 
        [p1 EulerDiff] = CalculatePolicy(p0);
        err = max(max(abs(p1 - p0))); %update e 
        disp([err count]);
        p0 = p1;
        count = count + 1;
        save(fname,'p1', 'nwgrid','xigrid','B','rho','mean_eta','err');
    end;
    toc;
    if count == maxiter;
        disp('maximum iteration reached');
    end;
      % Calculate consumption and savings policy
    c = InverseMargU(p1,mu);
    a = repmat(nwgrid,1,xisize) - c;
    figure('Name',['rho = ' num2str(rho) ', B = ', num2str(B)]);
    subplot(1,2,1);mesh(nwgrid,xigrid,c');
    xlabel('Networth');
    ylabel('\xi_t');
    title('consumption');
    subplot(1,2,2);mesh(nwgrid,xigrid,a');
    xlabel('Networth');
    ylabel('\xi_t');
    title('savings');

end;



