function [sL,sR]=exact_speeds_FULL(hL,hR,uL,uR,cL,cR,gravit)


tol=1.0E-12;
niter=100;

%Hcrit=(uR-uL) -2*(cL+cR); %depth positivity condition


%if (hL<=0 || hR<=0 || Hcrit>=0 )
    %wet bed case
    %evaluate the flow depth and velocity in the star region
    [hS,uS]=wetbed(hL,hR,uL,uR,cL,cR,gravit,niter,tol);
    if hS>=hL %left wave is a shock wave
        qL=sqrt((hS+hL)*hS/(2*hL^2));
        sL=uL-cL*qL;
    else %left wave is a rarefaction
        sL=uL-cL;
    end
    if hS>=hR %right wave is a shock
        qR=sqrt((hS+hR)*hS/(2*hR^2));
        sR=uR+cR*qR;
    else %right wave is a rarefaction
        sR=uR+cR;
    end
    
%end    
    
       
end




















