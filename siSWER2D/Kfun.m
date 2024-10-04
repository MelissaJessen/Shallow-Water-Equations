% Hydrolic conductivity given by the Genuchten constitutive relation
function K=Kfun(psi)
global alpha thetas thetar n m Ks psic 
sat = (Theta(psi)-thetar)/(thetas-thetar); % saturation 
K = Ks*sqrt(sat).*(1-(1-sat.^(1/m)).^m).^2;
