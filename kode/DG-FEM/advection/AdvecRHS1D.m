 function rhsu = AdvecRHS1D(u,param,maps,varargin)
% Strong form of advection equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = param.a;
K = param.K;
rx = param.rx;
Dr = param.Dr;
nx = param.nx;
Fscale = param.Fscale;
alpha = 1;
surfint = param.surfint;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;

jalpha = (1-0.5*alpha);
du = zeros(2,K);



for i = 1:K
    um  = u(vmapM(2*i-1:2*i));
    up  = u(vmapP(2*i-1:2*i));
    
    du(1,i) = -a*(um(1)-up(1))*jalpha;
    du(2,i) = a*(um(2)-up(2))*jalpha;
end

time = varargin{:};



if ~param.periodic
    uin = param.usol(0,time,a);
    du (mapI) = -a*(u(vmapI)-uin)*jalpha;
    du (mapO) = 0;
end



rhsu = -a * rx .* (Dr * u) + surfint * (Fscale .* du);


end