function c = InverseMargU(MargU,mu)
%==========================================================================
% Find the inverse of the marginal utility function
%==========================================================================
c = MargU.^(-1/mu);