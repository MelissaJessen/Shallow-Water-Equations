dmin = min(dl,dr);
ds = (1.0/g) * (0.5*(cl+cr) - 0.25*(ur-ul))^2;
% if ds<=dmin
% 
% else
if ds>dmin
% if ds>=dmin % first one, logically not right
    gel = sqrt(0.5*g*(ds+dl)/(ds*dl));
    ger = sqrt(0.5*g*(ds+dr)/(ds*dr));
    ds = (gel*dl + ger*dr - (ur-ul))/(gel+ger);
end
    
    