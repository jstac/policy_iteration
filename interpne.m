    function ZI= interpne(X,Y,Z,XI,YI,method)
    
    % Huiyu Li 2013/11/07
    % interp2 but extrapolate using the nearest value on the grid
    XI = max(min(X),XI);
    YI = max(min(Y),YI);
    XI = min(max(X),XI);
    YI = min(max(Y),YI);
    ZI= interp2(X, Y, Z, XI,YI,method); 
  