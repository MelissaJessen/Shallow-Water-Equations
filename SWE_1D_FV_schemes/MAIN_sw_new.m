% This code solves the Saint-Venant equations (1D problem)
% using different finite volume methods.
%
% Contacts:
% annunziato.siviglia@unitn.it
%
%-----------------------------------------------------------------------
clear 
close all
clc

SMALL = 1.E-12;
a=0;                      %left boundary
b=50;                     %right boundary
gate=(a+b)/2;                  %gate position for the dambreak
%gate = 20;
% Condition: the gate must be within the domain [a,b]
if gate>b
    'Error gate>b';
    stop
end



%-------------------------------------------------------------------
% Discretization
mcells=100;    %number of cells for numerical simulations

% FLUX CHOICE
%1--> Godunov method with exact Rieamnn problem
%2--> Lax-Friedrichs
%3--> Lax-Wendroff
%4--> FORCE
%5--> HLLC
%6--> Flux-splitting (2020) UPWIND: this is our method!!!
iflux=1;


gravit=9.8;       %gravity
ntmaxi=100000;     %maximum number of time cycles
%timeout=6;       %final time for computation
cfl=0.9;           %Courant number (choose cfl <1 for stable solutions)
% Initial data for simulations

%Initial conditions for the Riemann problem
%Choose TEST to run
% TEST 1:  Left rarefaction - Right shock
% TEST 2:  Left rarefaction - Right shock
% TEST 3:  Two-shocks
% TEST 4:  Two rarefactions
% TEST 5:  Two rarefactions (Red sea)

itest = 3;

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
elseif itest == 8
    % DATA GENERATION
    gate = 50.0;
    h_init_L=1.0;
    h_init_R=1.0;  
    UL =0.0;
    UR =0.0;
    timeout=4.0;
end



%----------------------------------------------------------------------

dx=(b-a)/mcells;          %spatial step
igate=round((gate-a)/dx,0);        %index (integer) at which the gate is positioned
x=a:dx:b; %interface coordinate
xc=x(1:mcells)+dx/2; %cell centre coordinate

%set the initial condition
h=zeros(1,mcells);
Q=zeros(1,mcells);
h(1:igate)=h_init_L;
Q(1:igate)=h_init_L*UL;
h(igate+1:mcells)=h_init_R;
Q(igate+1:mcells)=h_init_R*UR;
flux=zeros(2,mcells+1);
flux_ADV=zeros(2,mcells+1);
flux_PREX=zeros(2,mcells+1);

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

    %-------------Calculate numerical fluxes-------
    for i=2:mcells %internal interfaces
        hL  = h(i-1);
        hR  = h(i);
        uL=Q(i-1)/h(i-1);
        uR=Q(i)/h(i);
        qL=Q(i-1);
        qR=Q(i);
        aL=sqrt(gravit*hL);
        aR=sqrt(gravit*hR);

        if iflux==1
            %hstar and ustar are the solution of the Riemann problem at
            [hstar,ustar]=exrs_FULL(hL,hR,uL,uR,aL,aR,gravit);
            %Two-Rarefactions approximation
            %             hstar=(((aL+aR)/2+(uL-uR)/4)^2)/gravit;
            %             ustar = 0.5*(uL+uR) + aL-aR;
            %Computation of the Godunov flux
            flux(1,i)=hstar*ustar;
            flux(2,i)=ustar^2*hstar+1/2*gravit*hstar^2;

        elseif iflux==2  %Lax-Friedrichs
            FL(1)=hL*uL;
            FL(2)=uL^2*hL+0.5*gravit*hL^2;

            FR(1)=hR*uR;
            FR(2)=uR^2*hR+0.5*gravit*hR^2;


            flux(1,i) =0.5*(FL(1)+FR(1))-0.5*dx/dt*(hR-hL);
            flux(2,i) =0.5*(FL(2)+FR(2))-0.5*dx/dt*(hR*uR-hL*uL);


        elseif iflux==3   %LAX-Wendroff
            FL(1)=hL*uL;
            FL(2)=uL^2*hL+0.5*gravit*hL^2;

            FR(1)=hR*uR;
            FR(2)=uR^2*hR+0.5*gravit*hR^2;

            hLW = 0.5*(hL+hR)-0.5*dt/dx*(FR(1)-FL(1));
            QLW = 0.5*(uL*hL+uR*hR)-0.5*dt/dx*(FR(2)-FL(2));
            uLW = QLW/hLW;
            
            flux(1,i) = hLW*uLW;
            flux(2,i) =uLW^2*hLW+1/2*gravit*hLW^2;

        elseif iflux==4  %FORCE
            %LF
            FL(1)=hL*uL;
            FL(2)=uL^2*hL+0.5*gravit*hL^2;
           

            FR(1)=hR*uR;
            FR(2)=uR^2*hR+0.5*gravit*hR^2;
           

            FLF(1)=1/2*(FL(1)+FR(1))-1/2*dx/dt*(hR-hL);
            FLF(2)=1/2*(FL(2)+FR(2))-1/2*dx/dt*(hR*uR-hL*uL);


            %LW
            hLW = 0.5*(hL+hR)-0.5*dt/dx*(FR(1)-FL(1));
            QLW = 0.5*(uL*hL+uR*hR)-0.5*dt/dx*(FR(2)-FL(2));
            uLW = QLW/hLW;


            FLW(1) = hLW*uLW;
            FLW(2) = uLW^2*hLW+1/2*gravit*hLW^2;
            
            %FORCE
            flux(1,i) =0.5*(FLF(1)+FLW(1));
            flux(2,i) =0.5*(FLF(2)+FLW(2));



        elseif iflux==5   %HLL FLUX
            [SL,SR]=hll(hL,hR,uL,uR,aL,aR,gravit);
            ustar = 0.5*(uL+uR)+aL-aR; %two-rarefactions approx.

            FL(1)=hL*uL;
            FL(2)=uL^2*hL+0.5*gravit*hL^2;

            FR(1)=hR*uR;
            FR(2)=uR^2*hR+0.5*gravit*hR^2;

            Fhhl=(SR*FL-SL*FR+SR*SL*([hR,hR*uR]-[hL,hL*uL]))/(SR-SL);

            if SL>=0
                flux(1:2,i)=FL';
            elseif SR<=0
                flux(1:2,i)=FR';
            else
                flux(1:2,i)=Fhhl';
            end

        elseif iflux==6   %SPLITTING 
            
            FL(1)=hL*uL;
            FL(2)=uL^2*hL+1/2*gravit*hL^2;

            FR(1)=hR*uR;
            FR(2)=uR^2*hR+1/2*gravit*hR^2;

            %TWO-RAREFACTION APPROXIMATION
            [hstar,qstar]=two_rar_approx_PREX(hL,hR,qL,qR,gravit);

            flux_PREX(1) = qstar;
            flux_PREX(2) = 0.5*gravit*hstar^2;

            %Evaluation of the Advection flux

            if  qstar > 0
                flux_ADV(1) = 0;
                flux_ADV(2) = qstar*(uL);
            else
                flux_ADV(1) = 0;
                flux_ADV(2) = qstar*(uR);
            end
            %
            %
            %Computation of the flux
            F_SPLIT(1)=flux_ADV(1)+flux_PREX(1);
            F_SPLIT(2)=flux_ADV(2)+flux_PREX(2);

            flux(1:2,i)=F_SPLIT';

        end

    end

    %transmisive boundary conditions
    flux(:,1)=flux(:,2); %left boundary2
    flux(:,mcells+1)=flux(:,mcells); %right boundary

    %UPDATE FORMULA
    h  = h   - dt/dx*(flux(1,2:mcells+1)-flux(1,1:mcells));
    Q  = Q   - dt/dx*(flux(2,2:mcells+1)-flux(2,1:mcells));


    u= Q./h;
    Fr = u./sqrt(gravit*h);


    hmax=max(h)*1.2;
    hmin=min(h)*0.8;




    % %Plotting results
    % 
    % subplot(2,1,1)
    % plot(xc,h,'o')
    % subplot(2,1,2)
    % plot(xc,Q./h,'o')
    % drawnow

    time=time+dt;
    if time>=timeout
        %Solution of the exact Riemann problem to compare with the numerical
        % solution
        %
        mexact = 5000;
        %mexact = mcells-1;
        [hexact,uexact,dum,xexact]=exact_sol(mexact,h_init_L,h_init_R,UL,UR,0,0,gravit,timeout,a-gate,b-gate);
        %     sprintf('%0.16f',hexact(3947))
        %     sprintf('%0.16f',uexact(3947))

        umax=max(uexact)*1.2;
        umin=min(uexact);


        qexact = uexact.*hexact;
        qmax=max(Q)*1.2;
        qmin=min(Q);
        Frexact = uexact./sqrt(gravit.*hexact);
        Frmax=max(Frexact)*1.2;
        Frmin=min(Frexact);
       

        subplot(2,1,1)
        plot(xc,h,'o')
        hold on;
        plot(xexact+gate-dx/2,hexact,'k-','LineWidth',2)
        axis([a b hmin hmax]);
        title (['Water height at time = ', num2str(timeout),' s'])
        xlabel('x');
        ylabel('h');
        legend('Numerical solution', 'Exact solution');

        subplot(2,1,2)
        plot(xc,Q./h,'o')
        hold on;
        plot(xexact+gate-dx/2,uexact,'k-','LineWidth',2)
        axis([a b umin umax]);
        title (['Velocity at time = ', num2str(timeout),' s'])
        xlabel('x');
        ylabel('u');
        legend('Numerical solution', 'Exact solution');

        break
    end

end

% SAVE

% Save water height and velocity
u = Q./h;
xexact = xexact+gate-dx/2;

%folder = 'C:\Users\Matteo\Shallow-Water-Equations\data';
%filename = ['torotest', num2str(itest),'flux', num2str(iflux)];
%fullpath = fullfile(folder, filename);

%save(fullpath,'xc','h','u','xexact','hexact','uexact');


