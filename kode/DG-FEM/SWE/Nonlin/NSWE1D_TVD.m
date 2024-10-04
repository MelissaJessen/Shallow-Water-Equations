function [h,hu,hout,huout,tout] = NSWE1D_TVD(h, hu, x,FinalTime, factor,RHS,param,maps)



% Np = param.Np;
% K  = param.K;

g = param.g;
time = 0;

% Runge-Kutta residual storage  
% resu = zeros(Np,K); 
% reseta = zeros(Np,K);
% compute time step size
u = hu./h;

lambda = max(max(abs([u+sqrt(g*h/2) ; u-sqrt(g*h/2)])));
CFL = 0.2;
xmin = min(abs(x(1,:)-x(2,:)));
dt = CFL*xmin/lambda;
dt  = factor*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 



if nargout>2
hout = zeros(Nsteps,size(h,1),size(h,2));
hout(1,:,:) = h;
huout = zeros(Nsteps,size(hu,1),size(hu,2));
huout(1,:,:) = hu;
tout = zeros(1,Nsteps+1);

end

if param.SL
Limiter = param.Limiter;
h = Limiter(x,h,param); 
%h = Limiter(x,h + param.b,param)-param.b; 

hu= Limiter(x,hu,param); 

end


% outer time step loop 
for tstep=1:Nsteps
  % 3rd order SSP Runge-Kutta
  timelocal = time + dt;
  % SSP RK Stage 1.
  [rhsh,rhshu] = RHS(h,hu,param,maps,timelocal);

  h1  = h  + dt*rhsh;
  hu1 = hu + dt*rhshu;
  
  
  %   Limit fields
  if param.SL
    h1 = Limiter(x,h1,param); 
    %h = Limiter(x,h + param.b,param)-param.b; 

    hu1= Limiter(x,hu1,param);
  end
  
  [rhsh,rhshu] = RHS(h1,hu1,param,maps,timelocal);
  h2   = (3*h  + h1  + dt*rhsh )/4;
  hu2  = (3*hu + hu1 + dt*rhshu)/4;

%   Limit fields
  if param.SL
  h2  = Limiter(x,h2,param);
  hu2 = Limiter(x,hu2,param);
  end
  % SSP RK Stage 3.
  [rhsh,rhshu] = RHS(h2,hu2,param,maps,timelocal);
  h   = (h  + 2*h2  + 2*dt*rhsh )/3;
  hu  = (hu + 2*hu2 + 2*dt*rhshu)/3;

  % Limit solution
  if param.SL
  h  =Limiter(x,h,param);
  hu =Limiter(x,hu,param);
  end
  % Increment time
    time = time+dt;
    if nargout >2
    hout(tstep+1,:,:) = h;
    huout(tstep+1,:,:) = hu;
    tout(tstep+1) = time;
    end
        
end
return
