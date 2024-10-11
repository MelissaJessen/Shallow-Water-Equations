% DATA GENERATION
%
%
%-----------------------------------------------------------------------
clear 
close all
clc

a=0;                      %left boundary
b=50;                     %right boundary
gate=(a+b)/2;             %gate position for the dambreak
% Condition: the gate must be within the domain [a,b]
if gate>b
    'Error gate>b';
    stop
end

%Initial conditions for the Riemann problem
%Choose TEST to run
% TEST 1:  Left rarefaction - Right shock
% TEST 2:  Left rarefaction - Right shock
% TEST 3:  Two-shocks
% TEST 4:  Two rarefactions
% TEST 5:  Two rarefactions (Red sea)

itest = 1;

if itest==1
    % TORO TEST 1  
    gate = 10.0;
    h_init_L=1.0;
    h_init_R=0.1;
    UL =2.5;
    UR =0.0;
    timeout=7.0;
elseif itest == 2
    % TORO TEST 2  - Løses med flux 2, 4, 5
    gate = 25.0;
    h_init_L=1.0;
    h_init_R=1.0;
    UL =-5.0;
    UR =5.0;
    timeout=2.5;
 elseif itest == 3
    % TORO TEST 3 
    gate = 20.0;
    h_init_L=1.0;
    h_init_R=0.00005;  % add a little bit
    UL =0.0;
    UR =0.0;
    timeout=4.0;
 elseif itest == 4
    % TORO TEST 4 
    gate = 30.0;
    h_init_L=0.00005;  % add a little bit
    h_init_R=1.0;
    UL =0.0;
    UR =0.0;
    timeout=4.0;
 elseif itest == 5
    % TORO TEST 5  - Løses med flux 2, 4, 5
    gate = 25.0;
    h_init_L=0.1; 
    h_init_R=0.1;
    UL =-3.0;
    UR =3.0;
    timeout=5.0;
  elseif itest == 6
    % TORO TEST TRUE 
    gate = 20.0;
    h_init_L=1.0;
    h_init_R=0.0;  
    UL =0.0;
    UR =0.0;
    timeout=4.0;
  elseif itest == 7
    % TORO TEST TRUE 
    gate = 30.0;
    h_init_L=0.0;
    h_init_R=1.0;  
    UL =0.0;
    UR =0.0;
    timeout=4.0;
end

gravit=9.8;              %gravity
mcells = 100;
dx=(b-a)/mcells;          %spatial step
igate=round((gate-a)/dx,0);        %index (integer) at which the gate is positioned
x=a:dx:b; %interface coordinate
xc=x(1:mcells)+dx/2; %cell centre coordinate

mexact = 5000;
%set the initial condition
h=zeros(1,mcells);
Q=zeros(1,mcells);
h(1:igate)=h_init_L;
Q(1:igate)=h_init_L*UL;
h(igate+1:mcells)=h_init_R;
Q(igate+1:mcells)=h_init_R*UR;

ntmaxi=100000;     %maximum number of time cycles
cfl=0.9;           %Courant number (choose cfl <1 for stable solutions)

%Starting of the time marching procedure
time=0;
for t=1:ntmaxi %time marching precedure
    %-------------Courant condition----------------
    %Compute eigenvalues

    lambda(1,:)=Q./h-sqrt(gravit*h);
    lambda(2,:)=Q./h+sqrt(gravit*h);

    % %     %Compute maximum eigenvalue for the Courant condition
    lambdamax=max(max(abs(lambda)));
    % %     %set the time step dt
    dt=cfl*dx/lambdamax;
    if t<10
       dt=dt*0.2;
    end
    %-------------Courant condition----------------
    time=time+dt;
    if time>=timeout

        %Solution of the exact Riemann problem to compare with the numerical
        [hexact,uexact,~,xexact]=exact_sol(mexact,h_init_L,h_init_R,UL,UR,0,0,gravit,timeout,a-gate,b-gate);

        umax=max(uexact)*1.2;
        umin=min(uexact);

        qexact = uexact.*hexact;
        qmax=max(Q)*1.2;
        qmin=min(Q);
        Frexact = uexact./sqrt(gravit.*hexact);
        Frmax=max(Frexact)*1.2;
        Frmin=min(Frexact);
       
        break
    end

end

% Save water height and velocity
u = Q./h;
xexact = xexact+gate-dx/2;

folder = 'C:\Users\Matteo\Shallow-Water-Equations\dataFNO';
filename = ['test',num2str(itest),'_t=',num2str(timeout)];
fullpath = fullfile(folder, filename);

save(fullpath,'xexact','hexact','uexact');


