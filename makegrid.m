%==========================================================================
% Generates grid a with asize number of points and ranging from amin
% to amax. ascale=1 gives an equidistant grid. increasing ascale
% gridpoints in the lower values of grid. 
%==========================================================================

function a=makegrid(amin, amax, asize, ascale)

a=ones(asize,1);                          %initialise grid
adelta=(amax-amin)/((asize-1)^ascale);    %gap between grid points
for i=1:asize
    a(i)=amin+adelta*(i-1)^ascale;
end;    
 
