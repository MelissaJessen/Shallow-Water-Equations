function [rhsh,rhshu] = NSWE_RHS1D(h,hu,param,maps,timelocal)
g = param.g;
% K = param.K;
rx = param.rx;
Dr = param.Dr;
% nx = param.nx;
Fscale = param.Fscale;
surfint = param.surfint;
Flux_fun  = param.flux;
bx = param.bx;
% Minv = param.Minv;



u = hu./h;

F1 = hu;
F2 = hu.*u + 0.5 * g *h.^2;

[df1,df2] = Flux_fun(F1,F2,h,hu,param,maps,timelocal);


rhsh = -rx.*(Dr*F1) + surfint * (Fscale.*(df1)); %
rhshu = -rx.*(Dr*F2) + surfint * (Fscale.*(df2)) - (g*h.*bx);
 end