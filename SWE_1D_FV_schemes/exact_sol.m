function [hexact,uexact,coexact,x] = exact_sol(mcells,hL,hR,uL,uR,coL,coR,gravit,time,xmin,xmax)
                                  

%mcells  =5000;
tol     =1.E-16;
niter   =100;

h  = zeros(1,mcells+1);
u  = zeros(1,mcells+1);
co = zeros(1,mcells+1);

%hexact=zeros(1,mcells+1);
%uexact=zeros(1,mcells+1);

% compute celerity on left and right states
cL=sqrt(gravit*hL);
cR=sqrt(gravit*hR);

Hcrit=(uR-uL)-2*(cL+cR); %depth positivity condition
scrsz=get(0,'ScreenSize');
%f1=figure('Position',[scrsz(3)/4 10 3*scrsz(4)/4 scrsz(3)]);
if (hL<=0 || hR<=0 || Hcrit>=0 )
    %dry bed case
    if hL<=0
       shR=uR+cR;
       stR=uR-2*cR;
       dx=(xmax-xmin)/mcells;
       x=xmin:dx:xmax;
                 
           
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

           end
    else
        if hR<=0
           shL=uL-cL;
           stL=uL+2*cL;
           dx=(xmax-xmin)/mcells;
           x=xmin:dx:xmax;
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
        else
            shL=uL-cL;
            ssL=uL+2*cL;
            ssR=uR-2*cR;
            shR=uR+cR;
            dx=(xmax-xmin)/mcells;
            x=xmin:dx:xmax;
       
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
    
    dx=(xmax-xmin)/mcells;
    x=xmin:dx:xmax;
    
        S = x/time;
        for i=1:mcells+1
            if S(i)<uS
                co(i) = coL;
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
                co(i) = coR;
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
        
end        
for i=1:mcells+1
    hexact(i)=h(i);
    uexact(i)=u(i);
    coexact(i)=co(i);
    %x(i)=x(i)-25;
end

end

