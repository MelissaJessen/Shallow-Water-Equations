function [u,eta,uout,etaout,tout] = LinSWE1D(u, eta, x,FinalTime, factor,param,maps)



rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
     5161836677717.0/13612068292357.0 ...
     1720146321549.0/2090206949498.0  ...
     3134564353537.0/4481467310338.0  ...
     2277821191437.0/14882151754819.0];
 rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];
     
% Np = param.Np;
% K  = param.K;

time = 0;

% Runge-Kutta residual storage  
% resu = zeros(Np,K); 
% reseta = zeros(Np,K);
resu = 0;
reseta = 0;
% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.1; dt   = CFL/(2*pi)*xmin; dt = factor*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


if nargout>2
uout = zeros(Nsteps,size(u,1),size(u,2));
uout(1,:,:) = u;
etaout = zeros(Nsteps,size(eta,1),size(eta,2));
etaout(1,:,:) = eta;
tout = zeros(1,Nsteps+1);

end


% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        [rhsu,rhseta] = LinSWE_RHS1D(u,eta,param,maps,timelocal);
        
        resu = rk4a(INTRK)*resu + dt*rhsu;
        reseta = rk4a(INTRK)*reseta + dt*rhseta;

        u = u+rk4b(INTRK)*resu;
        eta = eta+rk4b(INTRK)*reseta;

    end
    % Increment time
    time = time+dt;
    if nargout >2
    uout(tstep+1,:,:) = u;
    etaout(tstep+1,:,:) = eta;
    tout(tstep+1) = time;
    end
        
end
return
