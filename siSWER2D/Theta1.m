function y = Theta1(psi)
global psic 


logic = (psi <= psic);
y = Theta(psi).*logic + (Theta(psic) + dTheta(psic)*(psi-psic)).*(~logic);

end
