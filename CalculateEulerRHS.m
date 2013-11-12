function MargUtilTmrw = CalculateEulerRHS(c,nw,xi0,p0)
% Given consumption level, net worth and xi today and policy function
% calculate expected marginal utility tomorrow  
% consumption, networth and xi0 can be vector of the same size
% solves for each triple 

global draw_eta draw_eps draw_u rho beta R nwgrid xigrid;
n = length(c);
MargUtilTmrw = zeros(n,1);
for j = 1:n; 
    xi1_vec = ((xi0(j)).^rho).*draw_u; %tmrw's xi, draw_u are draws from the law of motion for xi
    y1_vec = xi1_vec.*draw_eps + draw_eta; % tmrw's income 
    a1 = nw(j) - c(j); % tmrw's asset level
    networth1_vec = a1*R + y1_vec; %tmrw's net worth
    ptmrw = interpne(xigrid,nwgrid,p0,xi1_vec,networth1_vec,'linear'); % interpolate p0 to get tmrw's marginal utilities 
    MargUtilTmrw(j) = beta*R*mean(ptmrw); % RHS of Euler equation
                                                   % use diag because the
                                                   % interpolation is for a
                                                   % all permutations of
                                                   % xi1_vec and
                                                   % networth1_vec. Diag
                                                   % gives the values we
                                                   % want
end;