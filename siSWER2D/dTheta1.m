function y=dTheta1(psi)
global alpha thetas thetar n m Ks psic 

logic = (psi <= psic);
y = dTheta(psi).*logic + dTheta(psic)*(~logic);

end
% if(psi<=psic) 
%     % left of the critical value, take the original derivative 
%     y = dTheta(psi);  
% else
%     % on the right of the critical value, keep the maximum derivative 
%     y = dTheta(psic); 
% end