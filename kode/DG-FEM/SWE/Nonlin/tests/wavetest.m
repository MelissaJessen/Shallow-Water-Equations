function [u] = wavetest(x,a,z,delta,alpha,beta)
G = @(x,beta,z) exp(-beta*(x-z).^2);
F = @(x,alpha,a) sqrt(max(1-alpha^2*(x-a).^2));

u = zeros(size(x));
xmask1 = -0.8<=x & x<=-0.6;
xmask2 = -0.4<=x & x<=-0.2;
xmask3 = 0   <=x & x<= 0.2;
xmask4 = 0.4 <=x & x<= 0.6;

x1 = x(xmask1);
x2 = x(xmask2);
x3 = x(xmask3);
x4 = x(xmask4);

u(xmask1) = 1/6*(G(x1,beta,z-delta)+G(x1,beta,z+delta)+4*G(x1,beta,z));
u(xmask2) = 1;
u(xmask3) = 1-abs(10*(x3-0.1));                          %alpha? står | z i bogen
u(xmask4) = 1/6*(F(x4,alpha,a-delta)+G(x4,alpha,a+delta)+4*F(x4,alpha,a));



% u(-0.8<=x & x<=-0.6) = 1/6*(G(x,beta,z-delta)+G(x,beta,z+delta)+4*G(x,beta,z));
% u(-0.4<=x & x<=-0.2) = 1;
% u(0   <=x & x<= 0.2) = 1-abs(10*(x-0.1));                          %alpha? står | z i bogen
% u(0.4 <=x & x<= 0.6) = 1/6*(F(x,alpha,a-delta)+G(x,alpha,a+delta)+4*F(x,alpha,a));







