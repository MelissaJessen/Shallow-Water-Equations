function [hS,qS]=two_rar_approx_PREX(hL,hR,qL,qR,gravit)

tol = 1.0E-12;
cL = sqrt(gravit*hL);
cR = sqrt(gravit*hR);
%use two rarefaction solution as starting value
% cond_h_pos_1 = 0.75/sqrt(gravit)*(qL-qR);
% cond_h_pos_2 = 0.5*(hL^1.5+hR^1.5);
% diff = cond_h_pos_1 +cond_h_pos_2;
% 
% 
hS2 = (0.75/sqrt(gravit)*(qL-qR) + 0.5*(hL^1.5+hR^1.5))^2;
hS = hS2^(1./3.);
if  hS < tol
    hS = tol;
    %sprintf('%0.16f',uexact(3947))
end

            
qS = 0.5*(qL+qR) + sqrt(gravit)/3*(hL^1.5-hR^1.5);

if  abs(qS) < tol
    qS = sign(qS)*tol;
end

%qS = qR + fR;


end


