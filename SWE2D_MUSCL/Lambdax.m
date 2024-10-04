% characteristic velocities in the x-direction
function L = Lambdax(Q)
global g;

L = zeros(size(Q));
u = Q(2,:,:)./Q(1,:,:);
c = sqrt(g*Q(1,:,:));
L(1,:,:) = u - c;
L(2,:,:) = u;
L(3,:,:) = u + c;

end