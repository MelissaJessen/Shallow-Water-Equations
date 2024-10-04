function [d,u] = samrig(s,param)
% sample the solution for the case in which the right state is dry;
ul = param.ul;
cl = param.cl;
dl = param.dl;
dr = param.dr;
ur = param.ur;
g = param.g;
shl = ul-cl;

if s<=shl
    % sampling point lies to the left of the rarefaction
    d = dl;
    u = ul;
else
    stl = ul + 2*cl;
    if s<= stl
        % sampling point liesi nside the rarefaction
        u = (ul+2*cl+2*s)/3;
        c = (ul+2*cl-s)/3;
        d = c^2/g;
    else
        %sampling point lies in right drybed state
        d = dr;
        u = ur;
    end
end
end