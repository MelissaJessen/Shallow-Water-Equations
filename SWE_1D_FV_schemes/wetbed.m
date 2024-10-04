function [hS,uS]=wetbed(hL,hR,uL,uR,cL,cR,gravit,niter,tol)

%use two rarefaction solution as starting value
hS=(1/gravit)*(0.5*(cL+cR)-0.25*(uR-uL))^2;

h0=hS;
cha=1;

for i=1:niter
    [fL,fLd]=geofun(hS,hL,cL,gravit);
    [fR,fRd]=geofun(hS,hR,cR,gravit);
    hS = hS -(fL+fR+uR-uL)/(fLd+fRd);
    cha = abs(hS-h0)/(0.5*(hS+h0));
    if cha<tol
        break
    end
    if hS<0
        hS=tol;
    end
    h0=hS;
end

uS=.5*(uL+uR)+0.5*(fR-fL);

end


