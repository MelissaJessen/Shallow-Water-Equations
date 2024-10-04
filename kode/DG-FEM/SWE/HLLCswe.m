function [df1,df2] = HLLCswe(F1,F2,h,hu,param,maps,timelocal)
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
vmapL = zeros(2*K,1);
vmapR = zeros(2*K,1);
df1  = zeros(2,K);
df2  = zeros(2,K);
Fs1  = zeros(2,K);
Fs2  = zeros(2,K);
F1m = F1(vmapM); F2m = F2(vmapM);

vmapL(2:2:end) = vmapM(2:2:end);
vmapR(2:2:end) = vmapP(2:2:end);
vmapL(1:2:end-1) = vmapP(1:2:end-1);
vmapR(1:2:end-1) = vmapM(1:2:end-1);


SL = min([u(vmapL)-c(vmapL),u(vmapR)-c(vmapR)],[],2);
SR  = max([u(vmapL)+c(vmapL),u(vmapR)+c(vmapR)],[],2);

F1L = F1(vmapL); F2L = F2(vmapL); F1R = F1(vmapR); F2R = F2(vmapR);
hL  = h(vmapL);  hR  = h(vmapR);  huL = hu(vmapL); huR = hu(vmapR);
uL = u(vmapL); uR = u(vmapR);


Ss = (SL.*hR.*(uR-SR)-SR.*hL.*(uL-SL))./(hR.*(uR-SR)-hL.*(uL-SL));


hsL = hL.*(SL-uL)./(SL-Ss);
hsR = hR.*(SR-uR)./(SR-Ss);

husL = hsL.*Ss;
husR = hsR.*Ss;

% husL = hsL.*um; % Unsure which definition of flowrate is correct.
% husR = hsR.*up;

            
FsL1 = F1L + SL.*(hsL-hL);
FsL2 = F2L + SL.*(husL-huL);

FsR1 = F1R + SR.*(hsR-hR);
FsR2 = F2R + SR.*(husR-huR);


for i = 1:length(SL)
    if 0<=SL(i)
        Fs1(i) = F1L(i);
        Fs2(i) = F2L(i);
    elseif (SL(i)<= 0 && 0<=Ss(i))
        Fs1(i) = FsL1(i);
        Fs2(i) = FsL2(i);
    elseif (Ss(i)<= 0 && 0<=SR(i))
        Fs1(i) = FsR1(i);
        Fs2(i) = FsR2(i);
    else % SR<=0
        Fs1(i) = F1R(i);
        Fs2(i) = F2R(i);
    end
end

if param.BC
    % Left boundary: Boundary condition is the left state of Riemann

    % Specify boundary
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


    %     F1m     = F1(vmapI);
    %     F2m     = F2(vmapI);

    SL1 = min(uL1-cL1,uR1-cR1);
    SR1 = max(uL1 + cL1,uR1 + cR1);

    Ss1 = (SL1*hR1*(uR1-SR1)-SR1*hL1*(uL1-SL1))/(hR1*(uR1-SR1)-hL1*(uL1-SL1));

    %%%
    hsL1 = hL1*(SL1-uL1)/(SL1-Ss1);
    hsR1 = hR1*(SR1-uR1)/(SR1-Ss1);

    husL1 = hsL1*Ss1;
    husR1 = hsR1*Ss1;

    FsL11 = F1L1 + SL1*(hsL1-hL1);
    FsL21 = F2L1 + SL1*(husL1-huL1);
    FsR11 = F1R1 + SR1*(hsR1-hR1);
    FsR21 = F2R1 + SR1*(husR1-huR1);

    if 0<=SL1
        Fs1(mapI) = F1L1;
        Fs2(mapI) = F2L1;
    elseif (SL1<= 0 && 0<=Ss1)
        Fs1(mapI) = FsL11;
        Fs2(mapI) = FsL21;
    elseif (Ss1<= 0 && 0<=SR1)
        Fs1(mapI) = FsR11;
        Fs2(mapI) = FsR21;
    else % SR<=0
        Fs1(mapI) = F1R1;
        Fs2(mapI) = F2R1;
    end

    %% BC 2
    hR2  = param.h_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    huR2 = param.hu_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    uR2  = huR2/hR2;

    hL2  =  h(vmapO);
    huL2 = hu(vmapO);
    uL2  = u(vmapO);

%     
    F1L2 = F1(vmapO);
    F2L2 = F2(vmapO);
%     F1R2 = huL2;
%     F2R2 = huL2*uL2 + 0.5 * g *hL2^2;

    F1R2 = huR2;
    F2R2 = huR2*uR2 + 0.5 * g *hR2^2;

    cR2 = sqrt(g*hR2);
    cL2 = sqrt(g*hL2);


    SL2 = min(uL2-cL2,uR2-cR2);
    SR2 = max(uL2 + cL2,uR2 + cR2);

    Ss2 = (SL2*hR2*(uR2-SR2)-SR2*hL2*(uL2-SL2))/(hR2*(uR2-SR2)-hL2*(uL2-SL2));

    %%%
    hsL2 = hL2*(SL2-uL2)/(SL2-Ss2);
    hsR2 = hR2*(SR2-uR2)/(SR2-Ss2);

    husL2 = hsL2*Ss2;
    husR2 = hsR2*Ss2;
    
    FsL12 = F1L2 + SL2*(hsL2-hL2);
    FsL22 = F2L2 + SL2*(husL2-huL2);
    FsR12 = F1R2 + SR2*(hsR2-hR2);
    FsR22 = F2R2 + SR2*(husR2-huR2);

    if 0<=SL2
        Fs1(mapO) = F1L2;
        Fs2(mapO) = F2L2;
    elseif (SL2<= 0 && 0<=Ss2)
        Fs1(mapO) = FsL12;
        Fs2(mapO) = FsL22;
    elseif (Ss2<= 0 && 0<=SR2)
        Fs1(mapO) = FsR12;
        Fs2(mapO) = FsR22;
    else % SR<=0
        Fs1(mapO) = F1R2;
        Fs2(mapO) = F2R2;
    end


end


df1(:) = nx(:).*(F1m-Fs1(:));
df2(:) = nx(:).*(F2m-Fs2(:));

end