function [r,norm_r] = ResidualOuter(psi,bath,rhs)

    global Nx Nz dTheta12;

    r = zeros(Nz+1,Nx);
    dTheta12 = zeros(Nz+1,Nx);  % difference of dTheta1 - dTheta2
    Mpsi = MatVecProd_psi(psi);

    % compute the Outer residual, see eq. (30) in w1_d5.pdf
    r(1:Nz,:) = Theta(psi(1:Nz,:)      ) + Mpsi(1:Nz,:) - rhs(1:Nz,:);
    r(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1,:) - rhs(Nz+1,:);

    norm_r = sqrt( sum( sum(r.*r) ) );

end