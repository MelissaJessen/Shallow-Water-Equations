% characteristics of the SWE
function L = Lambda(q)
global Nx g;

u = q(2,:)./q(1,:);
h = q(1,:);

L = zeros(2,Nx);
L(1,:) = u - sqrt(g*h);
L(2,:) = u + sqrt(g*h);

end