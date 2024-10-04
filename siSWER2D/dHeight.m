% derivative of the nonlinear volume function 
% Input: free surface elevation (eta) and bathymetry (b)
% Output: derivative of the hieght function (dH) 
function dH = dHeight(eta,b) 

dH = 1.0*(eta - b>0);

end