% 
function Mpsi = MatVecProd_psi(psi)

    global dx dz dt Nx Nz Kx Kz dTheta12;
    
    fxm = zeros(Nz+1,Nx);
    fxp = zeros(Nz+1,Nx);
    fzm = zeros(Nz+1,Nx);
    fzp = zeros(Nz+1,Nx);

    % x fluxes
    fxm(:,1     ) = 0;  % no flux through the left border of the computational domain
    fxm(:,2:Nx  ) =-Kx(:,2:Nx).*( psi(:,2:Nx) - psi(:,1:Nx-1) )/dx;
    fxp(:,1:Nx-1) = fxm(:,2:Nx);
    fxp(:,Nx    ) = 0;  % no flux through the left border of the computational domain

    % z flux
    fzm(1     ,:) = 0;  % no flux through the left border of the computational domain
    fzm(2:Nz+1,:) =-Kz(2:Nz+1,:).*( psi(2:Nz+1,:) - psi(1:Nz,:) )/dz;
    fzp(1:Nz  ,:) = fzm(2:Nz+1,:);
    fzp(Nz+1  ,:) = 0;  % no flux through the left border of the computational domain

    % finite volume update
    Mpsi = dt/dx*( fxp - fxm ) + dt/dz*( fzp - fzm ) + dTheta12.*psi;

end
















