function [rhsh,rhshu] = NSWE_RHS1D_MMS(h,hu,param,maps,timelocal)
g = param.g;
rx = param.rx;
Dr = param.Dr;
Fscale = param.Fscale;
surfint = param.surfint;
Flux_fun  = param.flux;
bx = param.bx;
x = param.x;
L = param.L;
u = hu./h;

F1 = hu;
F2 = hu.*u + 0.5 * g *h.^2;

[df1,df2] = Flux_fun(F1,F2,h,hu,param,maps,timelocal);


rhsh = -rx.*(Dr*F1) + surfint * (Fscale.*(df1)) - exp(x/L)*sin(timelocal) + exp(x/L)*(sin(timelocal) + 3)/L;
rhshu = -rx.*(Dr*F2) + surfint * (Fscale.*(df2)) + (g*h.*bx) + exp(x/L)*cos(timelocal) + exp(x/L)*(sin(timelocal) + 3).^2/(L*(cos(timelocal) + 3)) + g*exp(x/L).^2*(cos(timelocal) + 3).^2/L;
end