function [avg,jump] = AvgJump(um,up,nx)
avg = 0.5*(um+up);
jump = nx(1)*um + nx(2)*up;
end