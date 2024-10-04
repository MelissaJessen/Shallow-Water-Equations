function [d,u,converged_flag] = wetbed(d,u,niter,param)
chalen = param.chalen;
mcells = param.mcells;
cl = param.cl;
cr = param.cr;
dr = param.dr;
dl = param.dl;
ul = param.ul;
ur = param.ur;

% fld,fld,fr,frd i dont think these need to be defined as they are computed
% each step in geofun;
gate = param.gate;
g = param.g;
timout = param.timout;
tol = param.tol;

starte
d0 = ds;
converged_flag = 0;
for i = 1:niter
    [fl,fld] = geofun(ds,dl,cl,param);
    [fr,frd] = geofun(ds,dr,cr,param);
    ds = ds-(fl+fr+ur-ul)/(fld+frd);
    cha = abs(ds-d0)/(0.5*(ds+d0));
    % missing soem write stuff but that shouldn't matter;
    
    if cha<=tol
        converged_flag = 1;
        break
    end
    
    if ds<0
        ds = tol;
    end
    d0 = ds;
end

if ~converged_flag
    fprint('not_converged')
    return
end
% If not above solution must have converged

us = 0.5*(ul + ur) + 0.5*(fr - fl);
cs = sqrt(g*ds);

for i = 1:mcells
    xcoord = i*chalen/mcells-gate;
    s = xcoord/timout;
    % Sample solution throughout the wave structure at time timout
    [dsam,usam] = samwet(s,ds,us,dl,ul,dr,ur,cl,cr,cs,param);
    d(i) = dsam;
    u(i) = usam;
end
    
        

    
    
end