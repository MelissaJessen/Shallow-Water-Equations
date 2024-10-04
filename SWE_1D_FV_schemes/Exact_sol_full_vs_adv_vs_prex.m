% This code solves the exact RP of the 1D:
%  1) shallow water system
%  2) advection system
%  3) pressure system
% 
% Contacts: 
% annunziato.siviglia@unitn.it
%
%-----------------------------------------------------------------------
clear all
close all
clc

SMALL = 1.E-12;
a=-10;                      %left boundary
b=10;                     %right boundary
gate=(a+b)/2;                  %gate position for the dambreak
% Condition: the gate must be within the domain [a,b]
if gate>b
    'Error gate>b'
    stop
end



itest=2;   %test choice:

                % TEST 1: RIGHT RAREFACTION - LEFT SHOCK
                % TEST 2: RIGHT SHOCK - LEFT SHOCK
                % TEST 3: RIGHT RAREFACTION - LEFT RAREFACTION (RED-SEA)
                % TEST 4: STATIONARY CONTACT

                %TESTS NOT INCLUDED IN THE PAPER
                % TEST 5:  Left rarefaction - Right shock
                % TEST 6:  Two rarefactions
                % TEST 7:  Isolated right shock (supersonic)
                % TEST 8:  Isolated left rarefaction 

if itest==1
    %TEST 1: RIGHT RAREFACTION - LEFT SHOCK
     gate=(a+b)/2;     
     h_init_L=1.0;
     h_init_R=0.6944444;
     UR=2.515045;
     UL=UR;
    coL_init = 1.0; %pollutant concentration
    coR_init = 0.0;
    timeout=0.5;  
    elseif itest == 2  
% TEST 2: RIGHT SHOCK - LEFT SHOCK (SLOWLY MOVING RIGHT SHOCK)
    h_init_L=.51;                    
    h_init_R=0.48;
    UL=2.5;
    UR=-5.8; %slowly-moving shock
    %UR=-2.3;
    coL_init = 1.0; %pollutant concentration
    coR_init = 0.0;
%     h_init_L=1;                    
%     h_init_R=0.5;
%     UL=100.0;
%     UR=0.0;
%     %UR=-5.3; %slowly-moving shock
    timeout=2;  
    elseif itest==3
   % TEST 3: RIGHT RAREFACTION - LEFT RAREFACTION (RED-SEA)
    h_init_L=1.0;                    
    h_init_R=1.0;
    UL =-3.0;
    UR = 3.0;
    coL_init = 1.0;
    coR_init = 0.0;
    timeout=1;  
    elseif itest==4
    %TEST 4: STATIONARY CONTACT
    h_init_L=1.0;                    
    h_init_R=1.0;
    UL = 0.0;
    UR = 0.0;
    coL_init = 1.0;
    coR_init = 0.55;
    timeout=6;  
elseif itest ==5
%     h_init_L=2.0;
%     h_init_R=0.1;
%     %h_init_R=0.0025;
%     UL=0.5;
%     UR=4.0;
%     coL_init = 1.0; %pollutant concentration
%     coR_init = 0.0;
%     timeout=0.5;  
    

%         a = -2;
%         b= 8;
%         gate=0;  
        
 %RIGHT WAVE
         a = -2;
         b= 8;
         gate=0;
         
%         h_init_L=1.0;
%         h_init_R=1.5;
%         UL=1.0;
%         UR=1.2;
%         h_init_L=1.0;
%         h_init_R=1.5;
%         UL=1.0;
%         UR=0.5;
%         h_init_L=1.0;
%         h_init_R=1.0;
%         UL=-0.5;
%         UR=-1.0;
%         h_init_L=1.0;
%         h_init_R=0.5;
%         UL=1.0;
%         UR=3.0;

%LEFT WAVE
        a = -3;
        b= 3;
        gate=0;
        
        
%         %Steady shock (UR<0 this cond. is from left shock)
%          h_init_L = 4.0;
%          h_init_R = 1.0;
%          UR = -2.0;
%          UL = -1.0;
         
         h_init_L = 3.5E-4;
         h_init_R = 0.3*h_init_L;
         UL = 1.0;
         UR = 2.0;
         %UL = -alpha*UR*sqrt(h_init_R/h_init_L);
        %Steady shock (UL>0 this cond. is from right shock)
%          h_init_L=1.0;
%          h_init_R=1.5;
%          UL = 1.0;
%          UR = -UL*sqrt(h_init_L/h_init_R);
         
         
         %UL = UR*sqrt(h_init_R/h_init_L);

        %Left shock
%          h_init_L=1.0;
%          h_init_R=1.5;
%          UL=-1.0;
%          UR=-1.1;
%Left rarefac
%          h_init_L=1.0;
%          h_init_R=0.5;
%          UL=-1.0;
%          UR=-1.1;
%Left shock
%          h_init_L=0.5;
%          h_init_R=1.0;
%          UL=-1.1;
%          UR=-1.0; 
%Left rarefac
%          h_init_L=1.0;
%          h_init_R=0.5;
%          UL=-1.5;
%          UR=-1.0;
        
        %Steady shock
%          h_init_L=1.0;
%          h_init_R=1.0;
%          UL = 1.0;
%          UR =-1.0;
        
         %Right shock
%          h_init_L=1.5;
%          h_init_R=1.0;
%          UL = 1.0;
%          UR =-1.0;
         
         %Left shock
%          h_init_L=0.5;
%          h_init_R=1.0;
%          UL = 1.0;
%          UR =-1.0;

%      Rarefaction centered at x=0
%          h_init_L=1.0;
%          h_init_R=1.0;
%          UL = -1.0;
%          UR = 1.0;
         
         %      Rarefaction centered at x=0
%          h_init_L=1.0;
%          h_init_R=1.0;
%          UL = -1.0;
%          UR = 1.0;

%             h_init_L=1.0;
%             h_init_R=1.0;
%             UL = -1.0;
%             UR = -1.0;
        
        coL_init = 1.0; %pollutant concentration
        coR_init = 0.0;
        timeout=0.05;  
        ustar_s = UL*sqrt(h_init_L/h_init_R);
        Speed = UR + ustar_s; 

%     a = -2;
%     b= 10;
%     gate=(a+b)/2;   
%     h_init_L=0.1;
%     h_init_R=0.1;
%     %h_init_R=0.0025;
%     UL=1.579404507382803;
%     UR=4.0;
%     coL_init = 1.0; %pollutant concentration
%     coR_init = 0.0;
%     timeout=0.5;  
  elseif itest ==6
    %TWO-RAREFACTIONS
    h_init_L=.5;                    
    h_init_R=0.51;
    UL=-0.3;
    UR=2.5;
    coL_init = 1.0; %pollutant concentration
    coR_init = 0.0;
    timeout=1.0;  
    elseif itest==7
        % Isolated shock 0.001
    gate=(a+b)/2;     
     h_init_L=0.066829783416185;
     h_init_R=0.001;
     UL=4.642433106032755;
     UR=0.0;
    coL_init = 1.0; %pollutant concentration
    coR_init = 0.0;
    timeout=2;  
    % Isolated shock 0.01
%     gate=(a+b)/2;     
%      h_init_L=0.171178918706455;
%      h_init_R=0.01;
%      UL=3.670582335681274;
%      UR=0.0;
%     coL_init = 1.0; %pollutant concentration
%     coR_init = 0.0;
%     timeout=2;
 elseif itest==8
     % Isolated left rarefaction 
     gate=(a+b)/2;
     h_init_L=1.0;
     h_init_R=0.171178918706455;
     UL=0.0;
     UR=3.670582335681274;
     coL_init = 1.0; %pollutant concentration
     coR_init = 0.0;
     timeout=2;

end


                   
                   
gravit=9.8;       %gravity 

   %Plotting results 

    
    
    time=timeout;
    
    %Solution of the exact Riemann problem to compare with the numerical
    % solution
    %
    mexact = 5000;
    %mexact = mcells-1;
    [hexact,uexact,coexact,xexact]=exact_sol(mexact,h_init_L,h_init_R,UL,UR,coL_init,coR_init,gravit,timeout,a-gate,b-gate);
    [hexact_adv,uexact_adv,coexact,xexact]=exact_sol_adv(mexact,h_init_L,h_init_R,UL,UR,coL_init,coR_init,gravit,timeout,a-gate,b-gate);
    %     sprintf('%0.16f',hexact(3947))
    %     sprintf('%0.16f',uexact(3947))
    
%     umax_ex=max(uexact);
%     umax=max(Q./h);
%     umax = max(umax,umax_ex);
%     
%     umin_ex=min(uexact);
%     umin=min(Q./h);
%     umin = min(umin,umin_ex);
%     
%     if (abs(umin-umax))<1.0E-9
%         umin = umin-0.1;
%         umax = umax +0.1;
%     end
%     
%     
     qexact = uexact.*hexact;
     qexact_adv = uexact_adv.*hexact_adv;
%     qmax = max(Q);
%     qmin = min(Q);
%     Frexact = uexact./sqrt(gravit.*hexact);
%     Frmax=max(Frexact)*1.2;
%     Frmin=min(Frexact);
%     
%     if (abs(qmin-qmax))<1.0E-9
%         qmin = qmin - 0.1;
%         qmax = qmax + 0.1;
%     end
%     
    hold off;
    subplot(3,1,1)
    plot(xexact+gate,hexact,'k-','LineWidth',2)
    hold on
    plot(xexact+gate,hexact_adv,'r-','LineWidth',2)
    %axis([a b]);
    %axis([a b hmin hmax]);
    title (['time = ', num2str(timeout),' s'])
    legend({'FULL','ADV'},'Location','southwest')
    xlabel('x');
    ylabel('h');
    
    hold off;
    subplot(3,1,2)
    plot(xexact+gate,uexact,'k-','LineWidth',2)
    hold on
    plot(xexact+gate,uexact_adv,'r-','LineWidth',2)
    %axis([a b]);
    %axis([a b umin umax]);
    title (['time = ', num2str(timeout),' s'])
    xlabel('x');
    ylabel('u');
    
    hold off;
    subplot(3,1,3)
    plot(xexact+gate,qexact,'k-','LineWidth',2)
    hold on
    plot(xexact+gate,qexact_adv,'r-','LineWidth',2)
    %axis([a b]);
    %axis([a b qmin qmax]);
    title (['time = ', num2str(timeout),' s'])
    xlabel('x');
    ylabel('q');
    
    
    

    

