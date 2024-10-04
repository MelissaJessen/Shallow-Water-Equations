% physical flux for the2D SWE in the x-direction
function F = FluxX(Q)
global g;

F = zeros(size(Q));

F(1,:,:) = Q(2,:,:);
F(2,:,:) = Q(2,:,:).^2./Q(1,:,:) + 0.5*g*Q(1,:,:).^2;
F(3,:,:) = Q(2,:,:).*Q(3,:,:)./Q(1,:,:);

end