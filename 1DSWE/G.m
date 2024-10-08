function res = G(hstar,hL,hR,uL,uR)
    res = Phi(hstar,hL) + Phi(hstar,hR) + uR - uL;
end