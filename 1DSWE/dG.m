function res = dG(hstar,hL,hR,uL,uR)
eps = 1e-8;

res = (G(hstar+eps,hL,hR,uL,uR) - G(hstar-eps,hL,hR,uL,uR))/(2*eps);

end