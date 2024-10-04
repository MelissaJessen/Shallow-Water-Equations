function [d,u] = sammid(s,param)
% sample the solution for the case in which the middle state is dry

ul = param.ul;
cl = param.cl;
ur = param.ur;
cr = param.cr;
dl = param.dl;
dr = param.dr;
g = param.g;

shl = ul - cl;
ssl = ul + 2*cl;
ssr = ur - 2*cr;
shr = ur + cr;

if s<=shl
    % sampling point lies the left of the left raefaction
    d = dl;
    u = ul;
% '1'
end

if (s>shl && s<=ssl)
    u = (ul + 2*cl + 2*s)/3;
    c = (ul + 2*cl - s)/3;
    d = c^2/g;
% '2'
end

if (s>ssl && s<=ssr)
    % sampling point lies in the middle dry bed region
    d = 0;
    u = 0;
% '3'
end

if (s>ssr && s<=shr)
    % sampling point lies inside the right rarefaction
    u = ( ur - 2*cr + 2*s)/3;
    c = (-ur + 2*cr + s)/3;
    d = c^2/g;
% '4'
end

if s>shr
    % sampling point to the right of the right rarefaction
    d = dr;
    u = ur;
% '5'
end

end


