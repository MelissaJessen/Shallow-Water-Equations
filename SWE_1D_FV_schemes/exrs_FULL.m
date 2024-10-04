function [hstar,ustar]=exrs_FULL(hL,hR,uL,uR,cL,cR,gravit)


tol=1.0E-8;
niter=100;

%Hcrit=(uR-uL) -2*(cL+cR); %depth positivity condition


%if (hL<=0 || hR<=0 || Hcrit>=0 )
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
    S = 0;
    if S<uS
       if left==1
          if S<sL
             hstar=hL;
             ustar=uL;
          else
             hstar=hS;
             ustar=uS;
          end
       else
          if S<shL
             hstar=hL;
             ustar=uL;
          else
             if S<stL
                ustar=(uL+2*cL+2*S)/3;
                C=(uL+2*cL-S)/3;
                hstar=C^2/gravit;
             else
                hstar=hS;
                ustar=uS;
             end
          end
       end
    else
       if right==1
          if S>sR
             hstar=hR;
             ustar=uR;
          else
             hstar=hS;
             ustar=uS;
          end
       else
          if S>shR
             hstar=hR;
             ustar=uR;
          else
             if S>stR
                ustar=(uR-2*cR+2*S)/3;
                C=(-uR+2*cR+S)/3;
                hstar=C^2/gravit;
             else
                hstar=hS;
                ustar=uS;
             end
          end
       end
    end
%end    
    
       
end




















