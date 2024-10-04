function [df1,df2] = HLLswe(F1,F2,h,hu,param,maps,timelocal)
g = param.g;
K = param.K;
% nx = param.nx;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;
% EToV  = param.EToV; 
% N = param.N;
nx = param.nx;
u   = hu./h;
c   = sqrt(g*h);
F1m = F1(vmapM); F2m = F2(vmapM);
vmapL = zeros(2*K,1);
vmapR = zeros(2*K,1);
df1  = zeros(2,K);
df2  = zeros(2,K);
Fs1  = zeros(2,K);
Fs2  = zeros(2,K);

vmapL(2:2:end) = vmapM(2:2:end);
vmapR(2:2:end) = vmapP(2:2:end);
vmapL(1:2:end-1) = vmapP(1:2:end-1);
vmapR(1:2:end-1) = vmapM(1:2:end-1);


SL = min([u(vmapL)-c(vmapL),u(vmapR)-c(vmapR)],[],2);
SR  = max([u(vmapL)+c(vmapL),u(vmapR)+c(vmapR)],[],2);

F1L = F1(vmapL); F2L = F2(vmapL); F1R = F1(vmapR); F2R = F2(vmapR);
hL  = h(vmapL);  hR  = h(vmapR);  huL = hu(vmapL); huR = hu(vmapR);


for i = 1:length(SL)
    if 0<=SL(i)
        Fs1(i) = F1L(i);
        Fs2(i) = F2L(i);
    elseif (SL(i)<= 0 && 0<=SR(i))
        Fs1(i) = (SR(i)*F1L(i) - SL(i)*F1R(i) + SL(i)*SR(i)*(hR(i)-hL(i)))/(SR(i)-SL(i));
        Fs2(i) = (SR(i)*F2L(i) - SL(i)*F2R(i) + SL(i)*SR(i)*(huR(i)-huL(i)))/(SR(i)-SL(i));
    else % SR<=0
        Fs1(i) = F1R(i);
        Fs2(i) = F2R(i);
    end
end



if param.BC
    hR1  = h(vmapI);
    huR1 = hu(vmapI);
    uR1  = u(vmapI);
    hL1  = param.h_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    huL1 = param.hu_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    uL1  = huL1/hL1;
    F1L1 = huL1;
    F2L1 = huL1*uL1 + 0.5 * g *hL1^2;
    F1R1 = F1(vmapI);
    F2R1 = F2(vmapI);

    cR1 = sqrt(g*hR1);
    cL1 = sqrt(g*hL1);

    SL1 = min(uL1-cL1,uR1-cR1);
    SR1 = max(uL1 + cL1,uR1 + cR1);

    if 0<=SL1
        Fs1(mapI) = F1L1;
        Fs2(mapI) = F2L1;
    elseif (SL1<= 0 && 0<=SR1)
        Fs1(mapI) = (SR1*F1L1 - SL1*F1R1 + SL1*SR1*(hR1-hL1))/(SR1-SL1);
        Fs2(mapI) = (SR1*F2L1 - SL1*F2R1 + SL1*SR1*(huR1-huL1))/(SR1-SL1);
    else % SR<=0
        Fs1(mapI) = F1R1;
        Fs2(mapI) = F2R1;
    end


%     Right boundary
    hL2   = h(vmapO);
    huL2  = hu(vmapO);
    uL2   = u(vmapO);
    hR2  = param.h_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    huR2 = param.hu_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    uR2 = huR2/hR2;
    F1R2   = huR2;
    F2R2   = huR2*uR2 + 0.5 * g *hR2^2;
    F1L2 = F1(vmapO);
    F2L2 = F2(vmapO);

    cR2 = sqrt(g*hR2);
    cL2 = sqrt(g*hL2);

    SL2 = min(uL2-cL2,uR2-cR2);
    SR2 = max(uL2 + cL2,uR2 + cR2);

    if 0<=SL2
        Fs1(mapO) = F1L2;
        Fs2(mapO) = F2L2;
    elseif (SL2<= 0 && 0<=SR2)
        Fs1(mapO) = (SR2*F1L2 - SL2*F1R2 + SL2*SR2*(hR2-hL2))/(SR2-SL2);
        Fs2(mapO) = (SR2*F2L2 - SL2*F2R2 + SL2*SR2*(huR2-huL2))/(SR2-SL2);
    else % SR<=0
        Fs1(mapO) = F1R2;
        Fs2(mapO) = F2R2;
    end


end

df1(:) = nx(:).*(F1m - Fs1(:));
df2(:) = nx(:).*(F2m - Fs2(:));
end