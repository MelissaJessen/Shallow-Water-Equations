function [d,u] = samlef(s,param)
% compute solution in the case when state is dry
ur = param.ur;
cr = param.cr;
g  = param.g;
dl = param.dl;
ul = param.ul;
dr = param.dr;
ur = param.ur;

shr = ur + cr;

if s>=shr
    % sampling point lies to the right of the rarefaction
    d = dr;
    u =  ur;
else
    str = ur-2*cr;
    if s>=str
        % sampling point lies inside the rarefaction
        u = (ur-2*cr + 2*s)/3;
        c = (-ur + 2*cr + s)/3;
        d = c^2/g;
    else
        %sampling point lies in drybed
        d = dl;
        u = ul;
    end
end