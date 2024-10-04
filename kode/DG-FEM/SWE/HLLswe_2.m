function [df1,df2] = HLLswe_2(F1,F2,h,hu,param,maps,timelocal)
% This version of the HLL flux is made to test whether we can abandon the
% left / right notation of the HLL flux in favor of solely relying on the
% +/- notation og the nodg framework
% 19/1-2022
g = param.g;
K = param.K;
% nx = param.nx;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;
nx = param.nx;
u   = hu./h;
c   = sqrt(g*h);
df1  = zeros(2,K);
df2  = zeros(2,K);
Fs1  = zeros(2,K);
Fs2  = zeros(2,K);


% 
% Sm = min([u(vmapM)-c(vmapM),u(vmapP)-c(vmapP)],[],2);
% Sp  = max([u(vmapM)+c(vmapM),u(vmapP)+c(vmapP)],[],2);


Sm = u(vmapM)-c(vmapM);
Sp  = u(vmapP)+c(vmapP);

F1m = F1(vmapM); F2m = F2(vmapM); F1p = F1(vmapP); F2p = F2(vmapP);
hm  = h(vmapM);  hp  = h(vmapP);  hum = hu(vmapM); hup = hu(vmapP);


% Maybe this can be done faster with masks?
for i = 1:length(Sm)
    if 0<=Sm(i)
        Fs1(i) = F1m(i);
        Fs2(i) = F2m(i);
    elseif (Sm(i)<= 0 && 0<=Sp(i))
%         Fs1(i) = nx(i) * (Sp(i)*F1m(i) - Sm(i)*F1p(i) + Sm(i)*Sp(i)*(hp(i)-hm(i)))/(Sp(i)-Sm(i));
%         Fs2(i) = nx(i) * (Sp(i)*F2m(i) - Sm(i)*F2p(i) + Sm(i)*Sp(i)*(hup(i)-hum(i)))/(Sp(i)-Sm(i));
        Fs1(i) = (Sp(i)*F1m(i) - Sm(i)*F1p(i) + nx(i) * Sm(i)*Sp(i)*(hp(i)-hm(i)))/(Sp(i)-Sm(i));
        Fs2(i) = (Sp(i)*F2m(i) - Sm(i)*F2p(i) + nx(i) * Sm(i)*Sp(i)*(hup(i)-hum(i)))/(Sp(i)-Sm(i));
        
%         Fs1(i) = (Sp(i)*F1m(i) - Sm(i)*F1p(i) + Sm(i)*Sp(i)*(hp(i)-hm(i)))/(Sp(i)-Sm(i));
%         Fs2(i) = (Sp(i)*F2m(i) - Sm(i)*F2p(i) + Sm(i)*Sp(i)*(hup(i)-hum(i)))/(Sp(i)-Sm(i));
    else % SR<=0
        Fs1(i) = F1p(i);
        Fs2(i) = F2p(i);
    end
end



if param.BC

    % Left boundary: Boundary condition is the left state of Riemann
    % problem

    hm1  = h(vmapI);
    hum1 = hu(vmapI);
    um1  = u(vmapI);
    hp1  = param.h_inflow(param.x(vmapI),timelocal,h,hu,param,maps); % Have to let this parameter be free otherwise we cant adjust inflow freely
    hup1 = param.hu_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    up1  = hup1/hp1;
    F1p1 = hup1;
    F2p1 = hup1*up1 + 0.5 * g *hp1^2;
    F1m1 = F1(vmapI);
    F2m1 = F2(vmapI);

    cm1 = sqrt(g*hm1);
    cp1 = sqrt(g*hp1);

    Sm1 = um1-cm1;
    Sp1 = up1 + cp1;

    if 0<=Sm1
        Fs1(mapI) = F1m1;
        Fs2(mapI) = F2m1;
    elseif (Sm1<= 0 && 0<=Sp1)               % At the very least there must be a sign switch here
%         Fs1(mapI) = nx(mapI) * (Sp1*F1m1 - Sm1*F1p1 + Sm1*Sp1*(hp1-hm1))/(Sp1-Sm1);
%         Fs2(mapI) = nx(mapI) * (Sp1*F2m1 - Sm1*F2p1 + Sm1*Sp1*(hup1-hum1))/(Sp1-Sm1);
        Fs1(mapI) = (Sp1*F1m1 - Sm1*F1p1 + nx(mapI) * Sm1*Sp1*(hp1-hm1))/(Sp1-Sm1);
        Fs2(mapI) = (Sp1*F2m1 - Sm1*F2p1 + nx(mapI) * Sm1*Sp1*(hup1-hum1))/(Sp1-Sm1);
%         Fs1(mapI) = (Sp1*F1m1 - Sm1*F1p1 + Sm1*Sp1*(hp1-hm1))/(Sp1-Sm1);
%         Fs2(mapI) = (Sp1*F2m1 - Sm1*F2p1 + Sm1*Sp1*(hup1-hum1))/(Sp1-Sm1);
    else % SR<=0
        Fs1(mapI) = F1p1;
        Fs2(mapI) = F2p1;
    end


%     Right boundary
    hm2  = h(vmapO);
    hum2 = hu(vmapO);
    um2  = u(vmapO);
    hp2  = param.h_inflow(param.x(vmapO),timelocal,h,hu,param,maps); % Have to let this parameter be free otherwise we cant adjust inflow freely
    hup2 = param.hu_inflow(param.x(vmapO),timelocal,h,hu,param,maps);
    up2  = hup2/hp2;
    F1p2 = hup2;
    F2p2 = hup2*up2 + 0.5 * g *hp2^2;
    F1m2 = F1(vmapO);
    F2m2 = F2(vmapO);

    cm2 = sqrt(g*hm2);
    cp2 = sqrt(g*hp2);

    Sm2 = um2-cm2;
    Sp2 = up2 + cp2;

    if 0<=Sm2
        Fs1(mapO) = F1m2;
        Fs2(mapO) = F2m2;
    elseif (Sm2<= 0 && 0<=Sp2)               % At the very least there must be a sign switch here
%         Fs1(mapO) = nx(mapO) * (Sp2*F1m2 - Sm2*F1p2 + Sm2*Sp2*(hp2-hm2))/(Sp2-Sm2);
%         Fs2(mapO) = nx(mapO) * (Sp2*F2m2 - Sm2*F2p2 + Sm2*Sp2*(hup2-hum2))/(Sp2-Sm2);
        Fs1(mapO) = (Sp2*F1m2 - Sm2*F1p2 + nx(mapO) * Sm2*Sp2*(hp2-hm2))/(Sp2-Sm2);
        Fs2(mapO) = (Sp2*F2m2 - Sm2*F2p2 + nx(mapO) * Sm2*Sp2*(hup2-hum2))/(Sp2-Sm2);
%         Fs1(mapO) = (Sp2*F1m2 - Sm2*F1p2 + Sm2*Sp2*(hp2-hm2))/(Sp2-Sm2);
%         Fs2(mapO) = (Sp2*F2m2 - Sm2*F2p2 + Sm2*Sp2*(hup2-hum2))/(Sp2-Sm2);
    else % SR<=0
        Fs1(mapO) = F1p2;
        Fs2(mapO) = F2p2;
    end

end

% df1(:) = nx(:).*(F1m - nx(:).*Fs1(:));
% df2(:) = nx(:).*(F2m - nx(:).*Fs2(:));

df1(:) = nx(:).*(F1m - Fs1(:));
df2(:) = nx(:).*(F2m - Fs2(:));
end

