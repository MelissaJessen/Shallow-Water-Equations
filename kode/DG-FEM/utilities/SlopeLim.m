function [ulim,xmid,aa,bb] = SlopeLim(x,u,param)
% vmapM = maps.vmapM;
% b     = param.b;
% M     = param.M;
Np    = param.Np;
K     = param.K;
V = param.V;
invV = param.invV;
rx = param.rx;
Dr = param.Dr;
J = param.J;
type = param.SLtype;



h = 2*J(1,:);
xmid = ones(Np,1)*(x(1,:) + h/2);

uh = invV*u; % Transfer to modal form


% Compute linear polynomial
ul = uh; ul(3:Np,:) = 0;ul = V*ul;
ulx = rx.*(Dr*ul);


% Compute cell averages
uh(2:Np,:)=0;
uavg = V*uh;

v = uavg(1,:);
vm = [v(1),v(1:K-1)];
vp = [v(2:K),v(K)]; 

aa = (v-vm)./h;
bb = (vp-v)./h;

psi = SlopeLimit(aa,bb,type,1,1);
ulim = uavg + (x-xmid).*psi;

% zlim = zhat;
end
