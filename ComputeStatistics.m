clear all; clc;

%==========================================================================
% Compares consumption behaviour from the income fluctuation model 
% when income is highly persistent versus when income has a permanent component
% 
% Also compares the effect of the tightness of the borrowing constraint
%
% Last updated 11/4/2013
%==========================================================================


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
% Simulate points for Monte Carlo integration
% Income process: 
% Y_t = xi_t * eps_t + eta_t
% eta_t ~ LN(mean_eta,sig_eta^2)
% ln xi_t ~ rho*ln(xi_{t-1}) + ln u_t, u_t ~LN (0,sig_xi^2)
% eps_t ~ LN(0,sig_eps^2)
sig_eps = 0.05;                 % std of multiplicative transitory shock
sig_u = 0.01;                  % std of persistent component
sig_xi0 = 0.15;                 % xi0 ~LN(0, sig_xi0^2)
sig_eta = 0.001;                % std of additive transitory shockmean_eta = 0; 

%Choose experiment
exercise = 1;                

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

fig = figure;

for rr=1:length(rho_vec);
    rho = rho_vec(rr);
    fname = ['Exercise_', num2str(exercise), '_rho',num2str(round(100*rho))];
    load(fname);    
    
    xisize = length(xigrid);
    nwsize = length(nwgrid);
    
    nwgridplot = min(nwgrid):1:25+min(nwgrid);
    xigridplot = 1:1:25;
    p1plot = interpne(xigrid,nwgrid,p1,xigridplot',nwgridplot,'linear');
   
    % Calculate consumption and savings policy
    c = InverseMargU(p1plot,mu);
    a = repmat(nwgridplot',1,length(xigridplot)) - c;

    nrow = 3;
    ncol = length(rho_vec);
    
    g1 = rr;
    g2 = rr + ncol;
    g3 = rr + ncol*2;

    figure(fig);
    subplot(nrow,ncol,g1); mesh(nwgridplot,xigridplot,c');
    xlabel('net worth','FontSize',12);
    ylabel('\xi_t','FontSize',12);
    zlabel('consumption','FontSize',12);
    axis([0 max(nwgridplot) 0 max(xigridplot) 0 30]);
    title(['\rho = ',num2str(rho)],'FontSize',18);
    subplot(nrow,ncol,g2); mesh(nwgridplot,xigridplot,a');
    xlabel('net worth','FontSize',12);
    ylabel('\xi_t','FontSize',12);
    zlabel('savings','FontSize',12);
    axis([0 max(nwgridplot) 0 max(xigridplot) 0 25]);
    
    % Simulate path of consumption   
    if exist('Panel');
        [N nVar Time] = size(Panel);
    else 
        N = 10000; Time = 2000;
        RandStream.setGlobalStream(s); %use the same random numbers for each parameter set
        c = InverseMargU(p1,mu); % consumption policy
        Panel = zeros(N,9,Time);
        a0 = 0;
        draws = randn(N,4);             % random draws for each shock
        draw_eta = exp(draws(:,1)*sig_eta); % simulated additive transitory shock  
        draw_eps = exp(draws(:,2)*sig_eps); % simulated multiplicative transitory shock  
        draw_u = exp(draws(:,3)*sig_u); % simulated xi_{t+1}/xi_t ratio
        draw_xi0 = exp(draws(:,4)*sig_xi0); %simulate initial draw of permanent shock


        for t = 1:Time        
            outgrid = sum(draw_xi0>max(xigrid));
            disp('Number of simulated xi out of grid');
            disp([t outgrid]);
            outgrid = sum(draw_xi0>max(nwgrid));
            disp('Number of simulated net worth out of grid');
            disp([t outgrid]);
            y =  draw_xi0.*draw_eps + draw_eta;
            nw0 =  a0*R + y;
            p = diag(interpne(xigrid,nwgrid,p1,draw_xi0,nw0','linear'));
            c = InverseMargU(p,mu);
            a = nw0 - c;
            bind = (a==-B);
            Panel(:,:,t) = [c nw0 a bind y draw_xi0 draw_eta draw_eps draw_u];
            a0 = a;
            draws = randn(N,4);             % random draws for each shock
            draw_eta = exp(draws(:,1)*sig_eta); % simulated additive transitory shock  
            draw_eps = exp(draws(:,2)*sig_eps); % simulated multiplicative transitory shock  
            draw_u = exp(draws(:,3)*sig_u); % simulated xi_{t+1}/xi_t ratio
            draw_xi0 = ((draw_xi0).^rho).* draw_u; % next period's xi

        end;
            
        save(fname,'p1', 'nwgrid','xigrid','B','rho','mean_eta','err','Panel');
    end;
    
    insur_coeff_mat = zeros(4,Time-1);

    for t = 2:Time;
        Change = log(Panel(:,:,t)) -  log(Panel(:,:,t-1));
        ChangeC = Change(:,1);
        temp =  cov([log(Panel(:,end-3:1:end,t)) ChangeC]) ; 
        insur_coeff = 1 - temp(:,end)./diag(temp); %permanent, transitory, transitory shock with change in consumption
        insur_coeff_mat(:,t-1) = insur_coeff(1:end-1);            
    end;

    figure(fig);
    subplot(nrow,ncol,g3); plot(2:1:Time, insur_coeff_mat(end-1:end,:));
    xlabel('period','FontSize',12);
    ylabel('insurability measure','FontSize',12);    
    legend('\epsilon_t','u_t','Location','SouthEast');
    axis([0,Time,-0.5,1]);
    clear Panel;
end;



