function [h,hu,uhout,uh2out,tout] = NSWE1D(h, hu, x,FinalTime, factor,RHS,param,maps)
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
g = param.g;
% xmin = min(abs(x(1,:)-x(2,:)));
% 
% CFL =0.9;
% dt  = CFL/(2*pi)*xmin;

% Compute time step size
u = hu./h;

lambda = max(max(abs([u+sqrt(g*h/2) ; u-sqrt(g*h/2)])));
CFL = 0.025; % For pipe case
% CFL = 1;
xmin = min(abs(x(1,:)-x(2,:)));
dt = CFL*xmin/lambda;
dt  = factor*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 


% if nargout>2
% uhout = zeros(Nsteps,size(h,1),size(h,2));
% uhout(1,:,:) = h;
% uh2out = zeros(Nsteps,size(hu,1),size(hu,2));
% uh2out(1,:,:) = hu;
% tout = zeros(1,Nsteps+1);
% 
% end

if nargout>2
uhout = zeros(size(h,1),size(h,2),Nsteps);
uhout(:,:,1) = h;
uh2out = zeros(size(hu,1),size(hu,2),Nsteps);
uh2out(:,:,1) = hu;
tout = zeros(Nsteps+1,1);

end

if param.SL
Limiter = param.Limiter;
end
% outer time step loop 
for tstep=1:Nsteps
    for INTRK = 1:5
        
        if param.SL
          h = Limiter(x,h,param);
%           h = Limiter(x,h + param.b,param) - param.b; 
          hu = Limiter(x,hu,param);
        end
        
        timelocal = time + rk4c(INTRK)*dt;
        [rhsh,rhshu] = RHS(h,hu,param,maps,timelocal);
        resu = rk4a(INTRK)*resu + dt*rhsh;
        reseta = rk4a(INTRK)*reseta + dt*rhshu;
        h = h+rk4b(INTRK)*resu;
        hu = hu+rk4b(INTRK)*reseta;
        

    end
    time = time+dt;

    if nargout >2
    uhout(:,:,tstep+1) = h;
    uh2out(:,:,tstep+1) = hu;
    tout(tstep+1) = time;
    end

        
end
return
