function [p1 EulerDiff] = CalculatePolicy(p0)

% given a guess of the price function, finds a new policy that satisfies
% the FOC

global pmin mu nwgrid xigrid cmax;

[n m] = size(p0);        % networth grid and shock grid sizes
pmax  = p0(1,1);         % marginal utility at the lowest asset and productivity level given by the initial guess
EulerDiff = zeros(n,m);  %difference between the LHS and RHS of the Euler equation at each grid point
p1    = pmin;            %initialize marginal utility at lowest marginal utility, that is initializing
                         %the consumption policy at the maximum feasible
                         %consumption level at each grid point. Borrowing
                         %constraint binds at these consumption levels

%for each asset level, compute the RHS of the Euler equation if the highest possible consumption
%level is chosen
MargUtilTmrwBnd = zeros(n,m);
for j = 1:m % for each xi today, calculate the RHS of Euler for each consumption and networth combination
    MargUtilTmrwBnd(:,j) = CalculateEulerRHS(cmax(:,j),nwgrid,repmat(xigrid(j),n,1),p0);
end;
[idx_n idx_xi] = find(MargUtilTmrwBnd > pmin); 
%find all states where the discounted marginal utility tomorrow is higher and so the 
%correct consumption level should be less than the maximum. 

%For states where the borrowing constraint does not bind, iterate to find the
%fixed point p1 = h(p1;p0), where h is the difference between the LHS and
%RHS of the Euler equation
numnbind = length(idx_n);  %number of states where the borrowing constraint does not bind      
for t = 1:numnbind         %for each nonbinding state
    i = idx_n(t); j = idx_xi(t); % index of networth and xi  
    nw = nwgrid(i);       %networth level today 
    xi0 = xigrid(j);      %xi today  
    
    % use secant method to find the solution
    phigh = pmax;         %Paper shows that the solution exists in this range
    plow  = pmin(i,j);  
    pguess1 = (phigh + plow)/2; %Secant uses two initial guess
    pguess0 = plow + (phigh-plow)/4;               
    
    % define the loss function
    % InverseMargU gives the consumption level for a price function
    loss = @(s) s - CalculateEulerRHS(InverseMargU(s,mu),nw,xi0,p0);
    
    loss1 = loss(pguess1);         
    loss0 = loss(pguess0);        

    maxiter = 5000; iter = 0;
    diff = 1;
    while diff > 1e-10&& iter < maxiter;
          iter = iter + 1;      
          update = loss1*(pguess1-pguess0)/(loss1-loss0);  % step for updating the guess
          pguess0 = pguess1; loss0 = loss1;
          pguess1 = pguess1 - update; 
          pguess1 = max(plow,pguess1);
          pguess1 = min(pmax, pguess1);
          loss1 = loss(pguess1);
          diff = abs(update);  % difference between guesses
    end;
    if iter >= maxiter; % signals maximum iteration reached 
       disp([diff i j]);   
    end;
    p1(i,j) = pguess1;
    EulerDiff(i,j) = loss1;
end;

