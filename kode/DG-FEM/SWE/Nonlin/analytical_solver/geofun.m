function [f,fd] = geofun(d,dk,ck,param)
g = param.g;
if d<=dk
    % rarefaction wave
    c = sqrt(g*d);
    f = 2.0*(c-ck);
    fd = g/c;
else
    ges = sqrt(0.5*g*(d+dk)/(d*dk));
    f   = (d-dk)*ges;
    fd = ges - 0.25*g*(d-dk) / (ges*d^2);
end
