% characteristic velocities in the y-direction
function L = Lambday(Q)
global g;

L = zeros(size(Q));
v = Q(3,:,:)./Q(1,:,:);
c = sqrt(g*Q(1,:,:));
L(1,:,:) = v - c;
L(2,:,:) = v;
L(3,:,:) = v + c;

end