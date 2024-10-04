function [df1,df2] = LaxFswe(F1,F2,h,hu,param,maps,timelocal)

g = param.g;
K = param.K;
nx = param.nx;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;

u = hu./h;
c = sqrt(g*h);


C = max(abs(u(vmapM)) + c(vmapM),abs(u(vmapP)) + c(vmapP));

df1  = zeros(2,K); 
df2  = zeros(2,K);


df1(:) = 0.5*nx(:).*(F1(vmapM) - F1(vmapP))-0.5*C.*(h(vmapM)-h(vmapP));
df2(:) = 0.5*nx(:).*(F2(vmapM) - F2(vmapP))-0.5*C.*(hu(vmapM)-hu(vmapP));

if param.BC
%     direct

    hm   = h(vmapI);
    hum  = hu(vmapI);
    hin  = param.h_inflow(param.x(vmapI),timelocal,h,hu,param,maps);
    huin = param.hu_inflow(param.x(vmapI),timelocal,h,hu,param,maps); 
%     huin = hum;
    uin = huin/hin;
    
    F1m    = F1(vmapI);
    F2m    = F2(vmapI);
    F1in   = hin*uin;
    F2in   = huin*uin + 0.5 * g *hin^2;
    
    cin = sqrt(g*hin);
    
    Cin = max(abs(u(vmapI)) + c(vmapI),abs(uin) + cin);

% 
    
    df1(mapI) = 0.5*nx(mapI)*(F1m-F1in) - 0.5*Cin*(hm-hin);
    df2(mapI) = 0.5*nx(mapI)*(F2m-F2in) - 0.5*Cin*(hum-huin);

% 
%%% right boundary direct
    
    hm      = h(vmapO);
    hum     = hu(vmapO);
    
    % Reflection??
%     hin  = hm;
%     huin = -hum;
%     uin  = huin/hin;  

    hin  = param.h_outflow(param.x(vmapO),timelocal,h,hu,param,maps);
    huin = param.hu_outflow(param.x(vmapO),timelocal,h,hu,param,maps); 
    uin  = huin/hin;  

    F1m   = F1(vmapO);
    F2m   = F2(vmapO);
    F1in  = hin*uin;
    F2in  = huin*uin + 0.5 * g *hin^2;
       
    cin = sqrt(g*hin);
    
    Cin = max(abs(u(vmapO)) + c(vmapO),abs(uin) + cin); % paper??
    
    df1(mapO) = 0.5*nx(mapO)*(F1m-F1in); - 0.5*Cin*(hm-hin);
    df2(mapO) = 0.5*nx(mapO)*(F2m-F2in); - 0.5*Cin*(hum-huin);    % new

end


end