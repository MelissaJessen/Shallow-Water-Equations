function [qk]=fun_qk(hs,hk)

y= hs/hk;
qk = sqrt(0.5*(y^2+y));

end
