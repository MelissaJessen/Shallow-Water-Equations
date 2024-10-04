% physical flux for the2D SWE in the x-direction
function G = FluxY(Q)
global g;

G = zeros(size(Q));

G(1,:,:) = Q(3,:,:);
G(2,:,:) = Q(2,:,:).*Q(3,:,:)./Q(1,:,:);
G(3,:,:) = Q(3,:,:).^2./Q(1,:,:) + 0.5*g*Q(1,:,:).^2;

end