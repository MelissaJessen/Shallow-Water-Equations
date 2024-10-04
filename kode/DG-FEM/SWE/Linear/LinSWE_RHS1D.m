function [rhsu,rhseta] = LinSWE_RHS1D(u,eta,param,maps,varargin)
% Strong form of advection equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g     = param.g;
A     = param.A;
h0    = param.h0;
omega = param.omega;
k     = param.k;


K = param.K;
rx = param.rx;
Dr = param.Dr;
nx = param.nx;
Fscale = param.Fscale;
% alpha = param.alpha;
surfint = param.surfint;
vmapM = maps.vmapM;
vmapP = maps.vmapP;
vmapI = maps.vmapI;
vmapO = maps.vmapO;
mapI  = maps.mapI;
mapO  = maps.mapO;
% R = param.R;
% Rinv = param.Rinv;
D = param.D;

du   = zeros(2,K);
deta = zeros(2,K);
du(:) = (u(vmapM)-u(vmapP));
deta(:) = (eta(vmapM)-eta(vmapP));

time = varargin{:};

C = abs(D(1));

% dfu   = 0.5*nx.*du*h0 - 0.5*C*(1-alpha)*nx.*(deta);
% dfeta = 0.5*nx.*deta*g - 0.5*C*(1-alpha)*nx.*(du);


dfu   = 0.5*nx.*du*h0 - 0.5*C*(deta);
dfeta = 0.5*nx.*deta*g - 0.5*C*(du);


if ~param.periodic
    % Left boundary
%     dirichlet:
    um        = u(vmapI);
    etam      = eta(vmapI);
    etain     = param.etasol(0,time,A,omega,k);
    uin       = param.usol(0,time,A,omega,k,h0);
    
    
    duin      = um-uin;
    detain    = etam-etain;
   
    
    dfu(mapI) = 0.5*nx(mapI)*h0*(duin) - 0.5*C*(detain);
    dfeta(mapI) = 0.5*nx(mapI)*g*(detain) - 0.5*C*(duin);
    
    
    % Right boundary
    um = u(vmapO);
    etam = eta(vmapO);
    
    % Dirichlet
%     etain  = param.etasol(param.xmax,time,A,omega,k);
%     uin    = param.usol(param.xmax,time,A,omega,k,h0);
    
%     uin = um;
%     etain = etam;
    uin = um;
    etain = -etam;

    
    % Transmission 
    duin   = um-uin;
    detain = etam-etain; 
    
    dfu(mapO) = 0.5*nx(mapO)*h0*(duin) - 0.5*C*(detain);

    dfeta(mapO) = 0.5*nx(mapO)*g*(detain) - 0.5*C*(duin);
    
%     % Reflection
%     dfu(mapO) = -C*etam;
%     dfeta(mapO) = nx(mapO)*etam;
%     

end


rhsu = -g*rx.*(Dr*eta) + surfint * (Fscale.*(dfeta));
rhseta = -h0*rx.*(Dr*u) + surfint * (Fscale.*(dfu));


end