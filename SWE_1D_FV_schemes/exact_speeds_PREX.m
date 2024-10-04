function [sL,sR]=exact_speeds_PREX(hL,hR,uL,uR,gravit)


tol=1.0E-6;
niter=50;

cL = sqrt(gravit*hL);
cR = sqrt(gravit*hR);
%Hcrit=(uR-uL) -2*(cL+cR); %depth positivity condition


%if (hL<=0 || hR<=0 || Hcrit>=0 )
    %wet bed case
    %evaluate the flow depth and velocity in the star region
    [hS,qS]=wetbed_PREX(hL,hR,uL,uR,cL,cR,gravit,niter,tol);
    
    if hS>=hL %left wave is a shock wave
        qL=sqrt(0.5*(hS/hL+1));
        sL=-cL*qL;
    else %left wave is a rarefaction
        sL = -cL;
    end
    if hS>=hR %right wave is a shock
        qR=sqrt(0.5*(hS/hR+1));
        sR=cR*qR;
    else %right wave is a rarefaction
        sR = cR;
    end

    
       
end




















