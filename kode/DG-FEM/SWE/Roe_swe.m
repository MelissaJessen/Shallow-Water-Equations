function [df1,df2] = Roe_swe(F1,F2,h,hu,param,maps,timelocal)
g = param.g;
K = param.K;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;
nx  = param.nx;
u   = hu./h;
c   = sqrt(g*h);
df1  = zeros(2,K);
df2  = zeros(2,K);
fs1  = zeros(2,K);
fs2  = zeros(2,K);



F1m = F1(vmapM); F2m = F2(vmapM); F1p = F1(vmapP); F2p = F2(vmapP);
hm  = h(vmapM);  hp  = h(vmapP);  
%hum = hu(vmapM); hup = hu(vmapP);


um = u(vmapM); up = u(vmapP); cm = c(vmapM); cp = c(vmapP);

% Roe averages
ubar = (um.*sqrt(hm) + up.*sqrt(hp))./(sqrt(hm) + sqrt(hp));
hbar = sqrt(hm.*hp);
cbar = sqrt(0.5*(cm.^2+cp.^2));

% Roe average eigenvalues / -vectors

lambda1bar = ubar-cbar;
lambda2bar = ubar + cbar;
K1 = ubar-cbar; K2 = ubar+cbar;
alpha1 = 0.5*(hp - hm - hbar./cbar.*(up-um));
alpha2 = 0.5*(hp - hm + hbar./cbar.*(up-um));



% -0.5*C.*(h(vmapM)-h(vmapP));
% -0.5*C.*(hu(vmapM)-hu(vmapP));


fs1(:) = 0.5*((F1m + F1p) - nx(:).*(alpha1 .* abs(lambda1bar) + alpha2 .* abs(lambda2bar)));
fs2(:) = 0.5*((F2m + F2p) - nx(:).*(alpha1 .* abs(lambda1bar) .* K1 + alpha2 .* abs(lambda2bar) .* K2));

df1(:) = 0.5*(nx(:).*(F1m - F1p) + (alpha1 .* abs(lambda1bar) + alpha2 .* abs(lambda2bar)));
df2(:) = 0.5*(nx(:).*(F2m - F2p) + (alpha1 .* abs(lambda1bar) .* K1 + alpha2 .* abs(lambda2bar) .* K2));



if param.BC
    % Left BC
    hm1  = h(vmapI);
%     hum1 = hu(vmapI);
    um1  = u(vmapI);
    hp1  = param.h_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    hup1 = param.hu_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    up1  = hup1/hp1;
    F1p1 = hup1;
    F2p1 = hup1*up1 + 0.5 * g *hp1^2;
    F1m1 = F1(vmapI);
    F2m1 = F2(vmapI);

    cm1 = c(vmapI);
    cp1 = sqrt(g*hp1);

    % Roe averages
    ubar1 = (um1*sqrt(hm1) + up1*sqrt(hp1))/(sqrt(hm1) + sqrt(hp1));
    hbar1 = sqrt(hm1*hp1);
    cbar1 = sqrt(0.5*(cm1^2+cp1^2));

    % Roe average eigenvalues / -vectors

    lambda1bar1 = ubar1-cbar1;
    lambda2bar1 = ubar1 + cbar1;
    K11 = ubar1-cbar1; K21 = ubar1+cbar1;
    alpha11 = 0.5*(hp1-hm1 - hbar1./cbar1*(up1-um1));
    alpha21 = 0.5*(hp1-hm1 + hbar1./cbar1*(up1-um1));
% 
%     fs1(mapI) = 0.5*((F1m1 + F1p1) - nx(mapI)*(alpha11 * abs(lambda1bar1) + alpha21 * abs(lambda2bar1)));
%     fs2(mapI) = 0.5*((F2m1 + F2p1) - nx(mapI)*(alpha11 * abs(lambda1bar1) * K11 + alpha21 * abs(lambda2bar1) .* K21));

    df1(mapI) = 0.5*(nx(mapI)*(F1m1 - F1p1) + (alpha11 * abs(lambda1bar1) + alpha21 * abs(lambda2bar1)));
    df2(mapI) = 0.5*(nx(mapI)*(F2m1 - F2p1) + (alpha11 * abs(lambda1bar1) * K11 + alpha21 * abs(lambda2bar1) .* K21));




%     Right boundary
    hm2  = h(vmapO);
%     hum2 = hu(vmapO);
    um2  = u(vmapO);
    hp2  = param.h_outflow(param.x(vmapO),timelocal,h,hu,param,maps); 
    hup2 = param.hu_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    up2  = hup2/hp2;
    F1p2 = hup2;
    F2p2 = hup2*up2 + 0.5 * g *hp2^2;
    F1m2 = F1(vmapO);
    F2m2 = F2(vmapO);

    cm2 = sqrt(g*hm2);
    cp2 = sqrt(g*hp2);

    
    % Roe averages
    ubar2 = (um2*sqrt(hm1) + up2*sqrt(hp1))/(sqrt(hm2) + sqrt(hp2));
    hbar2 = sqrt(hm2*hp2);
    cbar2 = sqrt(0.5*(cm2^2+cp2^2));

    % Roe average eigenvalues / -vectors

    lambda1bar2 = ubar2 - cbar2;
    lambda2bar2 = ubar2 + cbar2; % 
    K12 = ubar2-cbar2; K22 = ubar2+cbar2;
    alpha12 = 0.5*(hp2-hm2 - hbar2./cbar2*(up2-um2));
    alpha22 = 0.5*(hp2-hm2 + hbar2./cbar2*(up2-um2));

    df1(mapO) = 0.5*(nx(mapO) * (F1m2 - F1p2) + (alpha12 * abs(lambda1bar2) + alpha22 * abs(lambda2bar2)));
    df2(mapO) = 0.5*(nx(mapO) * (F2m2 - F2p2) + (alpha12 * abs(lambda1bar2) * K12 + alpha22 * abs(lambda2bar2) .* K22));



end

% df1(:) = nx(:).*(F1m - fs1(:));
% df2(:) = nx(:).*(F2m - fs2(:));
end



