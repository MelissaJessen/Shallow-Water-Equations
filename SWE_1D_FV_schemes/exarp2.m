%This code solves the Exact Riemann problem for the Saint-Venant
%equations. This has been developed at Univerity of Trento 
%by prof. E:F. Toro and Dr. G. Vignolia nd it is used in the Master course 
%of Numerical Methods for the Environment"
%Thi scode has been also been adopted in the course of Numerical modelling
%in hydraulics and river engineering at ETH. The tutor of this course is
%Dr. A. Siviglia.

clear all
close all
clc 
chalen  =50;
gate    =25;
gravit  =9.81;
mcells  =5000;
ntmaxi  =500000;
dt      =0.05;
tol     =1.E-6;
niter   =50;
timeout =2;
hL      =0.25;
uL      =0.0 ;
hR      =1.;
uR      =0.0;


h=zeros(1,mcells+1);
u=zeros(1,mcells+1);



% compute celerity on left and right states
cL=sqrt(gravit*hL);
cR=sqrt(gravit*hR);

Hcrit=(uR-uL)-2*(cL+cR); %depth positivity condition
scrsz=get(0,'ScreenSize');
f1=figure('Position',[scrsz(3)/4 10 3*scrsz(4)/4 scrsz(3)]);
if (hL<=0 || hR<=0 || Hcrit>=0 )
    %dry bed case
    if hL<=0
       shR=uR+cR;
       stR=uR-2*cR;
       xmin=1.1*(stR*timeout);
       xmax=1.1*(shR*timeout);
       dx=(xmax-xmin)/mcells;
       x=xmin:dx:xmax;
       for t=1:ntmaxi
           time=t*dt;
           S = x/time;
           for i=1:mcells+1
               if S(i) >= shR
                  h(i)=hR;
                  u(i)=uR;
               else
                  if S(i)>=stR
                     u(i)=(uR-2*cR+2*S(i))/3;
                     C =(-uR+2*cR+S(i))/3;
                     h(i)=C^2/gravit;
                  else
                     h(i)=hL;
                     u(i)=uL;
                  end
               end
%                if(timeout<=0)
%                    
%                    if(x(i)<gate)
%                        h(i)=hL;
%                        u(i)=uL;
%                    else
%                        h(i)=hR;
%                        u(i)=uR;
%                    end
%                end
           end
           subplot(3,1,1)
           plot(x+gate,h)
           xlabel('x')
           ylabel('h - flow depth')
           axis([xmin+gate xmax+gate, 0 1.1*max(hL,hR)]);
           subplot(3,1,2)
           plot(x+gate,u,'-g')
           xlabel('x')
           ylabel('u - flow velocity')
           axis([xmin+gate xmax+gate, 1.1*min(u) 1.1*max(u)]);
           subplot(3,1,3)
           for lc=1:11 %draw characteristic lines
               tc=max(0,x/(stR+(lc-1)/10*(shR-stR)) );
               plot(x+gate,tc,'-k')
               axis([xmin+gate xmax+gate, 0 1.1*timeout]);
               xlabel('x')
               ylabel('time')
               hold on
           end
           plot(x+gate,time*ones(1,mcells+1),'-r')
           hold off
           drawnow
           if time>=timeout
               break
           end
       end
    else
        if hR<=0
           shL=uL-cL;
           stL=uL+2*cL;
           xmin=1.1*(shL*timeout);
           xmax=1.1*(stL*timeout);
           dx=(xmax-xmin)/mcells;
           x=xmin:dx:xmax;
           for t=1:ntmaxi
               time=t*dt;
               S = x/time;
               for i=1:mcells+1
                   if S(i)<=shL
                       h(i)=hL;
                       u(i)=uL;
                   else
                       if S(i)<=stL
                           u(i)=(uL+2*cL+2*S(i))/3;
                           C=(uL+2*cL-S(i))/3;
                           h(i)=C^2/gravit;
                       else
                           h(i)=hR;
                           u(i)=uR;
                       end
                   end
               end
               subplot(3,1,1)
               plot(x+gate,h)
               xlabel('x')
               ylabel('h - flow depth')
               axis([xmin+gate xmax+gate, 0 1.1*max(hL,hR)]);
               subplot(3,1,2)
               plot(x+gate,u,'-g')
               xlabel('x')
               ylabel('u - flow velocity')
               axis([xmin+gate xmax+gate, 1.1*min(u) 1.1*max(u)]);
               axis([xmin+gate xmax+gate, 1.1*min(u) 1.1*max(u)]);
               subplot(3,1,3)
               for lc=1:11 %draw characteristic lines
                   tc=max(0,x/(shL+(lc-1)/10*(stL-shL)) );
                   plot(x+gate,tc,'-k')
                   xlabel('x')
                   ylabel('time')
                   axis([xmin+gate xmax+gate, 0 1.1*timeout]);
                   hold on;
               end
               plot(x+gate,time*ones(1,mcells+1),'-r')
               hold off
               drawnow
               if time>=timeout
                  break
               end
           end
        else
            shL=uL-cL;
            ssL=uL+2*cL;
            ssR=uR-2*cR;
            shR=uR+cR;
            xmin=1.1*(shL*timeout);
            xmax=1.1*(shR*timeout);
            dx=(xmax-xmin)/mcells;
            x=xmin:dx:xmax;
            for t=1:ntmaxi
                time=t*dt;
                S = x/time;
                for i=1:mcells+1
                    if S(i)<=shL
                        h(i)=hL;
                        u(i)=uL;
                    end
                    if S(i)>shL && S(i)<=ssL
                        u(i)=(uL+2*cL+2*S(i))/3;
                        C = (uL+2*cL-S(i))/3;
                        h(i)=C^2/gravit;
                    end
                    if S(i)>ssL && S(i)<=ssR
                        h(i)=0;
                        u(i)=0;
                    end
                    if S(i)>ssR && S(i)<=shR
                        u(i)=(uR-2*cR+2*S(i))/3;
                        C = (-uR+2*cR+S(i))/3;
                        h(i)=C^2/gravit;
                    end
                    if S(i)>shR
                        h(i)=hR;
                        u(i)=uR;
                    end
                end
                subplot(3,1,1)
                plot(x+gate,h)
                xlabel('x')
                ylabel('h - flow depth')
                axis([xmin+gate xmax+gate, 0 1.1*max(hL,hR)]);
                subplot(3,1,2)
                plot(x+gate,u,'-g')
                xlabel('x')
                ylabel('u - flow velocity')
                axis([xmin+gate xmax+gate, 1.1*min(u) 1.1*max(u)]);
                axis([xmin+gate xmax+gate, 1.1*min(u) 1.1*max(u)]);
                subplot(3,1,3)
                for lc=1:11 %draw characteristic lines
                    tc=max(0,x/(shL+(lc-1)/10*(ssL-shL)) );
                    tc1=max(0,x/(ssR+(lc-1)/10*(shR-ssR)) );
                    plot(x+gate,tc,'-k')
                    hold on
                    plot(x+gate,tc1,'-k')
                    xlabel('x')
                    ylabel('time')
                    axis([xmin+gate xmax+gate, 0 1.1*timeout]);
                end
                plot(x+gate,time*ones(1,mcells+1),'-r')
                hold off
                drawnow 
                if time>=timeout
                    break
                end
            end
        end
    end
else
    %wet bed case
    %evaluate the flow depth and velocity in the star region
    [hS,uS]=wetbed(hL,hR,uL,uR,cL,cR,gravit,niter,tol);
    cS=sqrt(gravit*hS);
    if hS>=hL %left wave is a shock wave
        qL=sqrt((hS+hL)*hS/(2*hL^2));
        sL=uL-cL*qL;
        left=1;
    else %left wave is a rarefaction
        shL=uL-cL;
        stL=uS-cS;
        left=0;
    end
    if hS>=hR %right wave is a shock
        qR=sqrt((hS+hR)*hS/(2*hR^2));
        sR=uR+cR*qR;
        right=1;
    else %right wave is a rarefaction
        shR=uR+cR;
        stR=uS+cS;
        right=0;
    end
    if left==1
        xmin=min(-5,1.1*sL*timeout);
    else
        xmin=min(-5,1.1*shL*timeout);
    end
    if right==1
        xmax=max(5,1.1*sR*timeout);
    else
        xmax=max(5,1.1*shR*timeout);
    end
%xmin=-15;
%xmax=15;
    dx=(xmax-xmin)/mcells;
    x=xmin:dx:xmax;
    for t=1:ntmaxi
        time=t*dt;
        S = x/time;
        for i=1:mcells+1
            if S(i)<uS
                if left==1
                    if S(i)<sL
                        h(i)=hL;
                        u(i)=uL;
                    else
                        h(i)=hS;
                        u(i)=uS;
                    end
                else
                    if S(i)<shL
                        h(i)=hL;
                        u(i)=uL;
                    else
                        if S(i)<stL
                            u(i)=(uL+2*cL+2*S(i))/3;
                            C=(uL+2*cL-S(i))/3;
                            h(i)=C^2/gravit;
                        else
                            h(i)=hS;
                            u(i)=uS;
                        end
                    end
                end
            else
                if right==1
                    if S(i)>sR
                        h(i)=hR;
                        u(i)=uR;
                    else
                        h(i)=hS;
                        u(i)=uS;
                    end
                else
                    if S(i)>shR
                       h(i)=hR;
                       u(i)=uR;
                    else
                        if S(i)>stR
                            u(i)=(uR-2*cR+2*S(i))/3;
                            C=(-uR+2*cR+S(i))/3;
                            h(i)=C^2/gravit;
                        else
                            h(i)=hS;
                            u(i)=uS;
                        end
                    end
                end
            end
        end
        subplot(3,1,1)
        plot(x+gate,h)
        xlabel('x')
        ylabel('h - flow depth')
        axis([xmin+gate xmax+gate 0 1.1*max(h)])
        subplot(3,1,2)
        plot(x+gate,u,'-g')
        xlabel('x')
        ylabel('u - flow velocity')
        axis([xmin+gate xmax+gate 1.1*min(u) 1.1*max(u)])
        
        
        max_time= 1.1*timeout;
        %max_time=3.5;
        subplot(3,1,3)
        if left==1
           tc=max(-0.1,x/sL);
           plot(x+gate,tc,'-b','LineWidth',3)
           axis([xmin+gate xmax+gate 0 max_time])
           hold on;
        else
           for lc=1:11
               tc=max(0,x/(shL+(lc-1)/10*(stL-shL)) );
               plot(x+gate,tc,'-k')
               axis([xmin+gate xmax+gate 0 max_time])
               hold on;
           end
        end
        if right==1
           tc=max(-0.1,x/sR);
           plot(x+gate,tc,'-b','LineWidth',3)
           axis([xmin+gate xmax+gate 0 max_time])
           hold on
        else
           for lc=1:11
               tc=max(0,x/(stR+(lc-1)/10*(shR-stR)) );
               plot(x+gate,tc,'-k')
               axis([xmin+gate xmax+gate 0 max_time])
               hold on
           end
        end
        plot(x+gate,time*ones(1,mcells+1),'-r')
        axis([xmin+gate xmax+gate 0. max_time])
        xlabel('x')
        ylabel('time')
        hold off;
        drawnow
        if time>=timeout
            break
        end
    end        
end




















