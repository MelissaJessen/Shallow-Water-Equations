function y=dTheta(psi)
global alpha thetas thetar n m Ks psic 


% unsaturated medium 
logic = (psi <= 0);
y = (alpha*n*m*(thetas - thetar)./((1+abs(alpha*psi).^n).^(m+1)).*abs(alpha*psi).^(n-1)).*logic; 
% saturated medium y = 0; 

end

% if(psi<=0)
%     % unsaturated medium 
%     y = alpha*n*m*(thetas - thetar)./((1+abs(alpha*psi).^n).^(m+1)).*abs(alpha*psi).^(n-1); 
% else
%     % saturated medium 
%     y = 0; 
% end