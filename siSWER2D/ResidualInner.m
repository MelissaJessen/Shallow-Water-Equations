function [r,norm_r] = ResidualInner(psi,psi_alpha,bath,rhs)

    global Nx Nz dTheta12;

    r = zeros(Nz+1,Nx);
    dTheta12 = zeros(Nz+1,Nx);  % difference of dTheta1 - dTheta2
    Mpsi = MatVecProd_psi(psi);

    % compute the Inner residual, see eq.(37) in w1_d5.pdf
    % inside of the orous medium:
    r(1:Nz,:) = Theta1(psi(1:Nz,:)) - Theta2(psi_alpha(1:Nz,:)) ...
               - dTheta2(psi_alpha(1:Nz,:)).*( psi(1:Nz,:) - psi_alpha(1:Nz,:) ) ...
               + Mpsi(1:Nz,:) - rhs(1:Nz,:);
    % inside of the Shallow water layer:
    r(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1,:) - rhs(Nz+1,:);

    norm_r = sqrt( sum( sum(r.*r) ) );

    dTheta12(1:Nz,:) = dTheta1(psi(1:Nz,:)) - dTheta2(psi_alpha(1:Nz,:));
    dTheta12(Nz+1,:) = dHeight(psi(Nz+1,:),bath);

end