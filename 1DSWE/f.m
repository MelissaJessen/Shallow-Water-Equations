% physical flux for the SWE
function F = f(q)

global g;

F = zeros(2,1);
F(1) = q(2);                          % h*u
F(2) = q(2)^2/q(1) + 0.5*g*q(1)^2;    % h*u^2 + 1/2*g*h^2

end