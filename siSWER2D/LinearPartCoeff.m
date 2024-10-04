% Compute the coefficients of the linear part Kx, Kz and the rhs
function [theta,rhs] = LinearPartCoeff(psi,theta,bath,ustar)
    global K Kx Kz Nx Nz dx dz dt g gamma Hx
    rhs =zeros(Nz+1,Nx);

    % compute volume fraction and hydraulic conductivity at the cell centers
    theta(1:Nz,:) = Theta( psi(1:Nz,:) );           % nonlinear volume fraction function inside of the porous medium
    theta(Nz+1,:) = Height( psi(Nz+1,:),bath );     % and in the shallow water layer
    K(1:Nz,:) = Kfun(psi(1:Nz,:));                  % Hydraulic cond inside of the porous medium
    K(Nz+1,:) = Kfun(psi(Nz+1,:));                  % Hydraulic cond inside of the shallow water layer

    % compute H at the x-faces in the shallow water layer
    H = theta(Nz+1,:);  % H at the cell centers
    Hx(1)    = 0;  % only suited for dry cells at the left boundary
    Hx(2:Nx) = max(0,max(H(2:Nx),H(1:Nx-1))); 
    Hx(Nx+1) = 0;  % only suited for dry cells at the left boundary
    Hx_tilde = Hx.^2./( Hx + dt*gamma + 1e-12);

    % Hydraulic conductivity at the x-faces, Kx
    % inside of the porous medium
    Kx(1:Nz,1)    = 0;
    Kx(1:Nz,2:Nx) = max(K(1:Nz,2:Nx),K(1:Nz,1:Nx-1));
    Kx(1:Nz,Nx+1) = 0;
    % inside of the shallow water layer
    Kx(Nz+1,:) = g*dt*Hx_tilde;

    % Hydraulic conductivity at z-faces, Kz
    % inside of the porous medium:
    Kz(1     ,1:Nx) = 0;
    Kz(2:Nz+1,1:Nx) = max(K(2:Nz+1,1:Nx),K(1:Nz,1:Nx));
    % inside of the shallow water layer
    Kz(Nz+2  ,1:Nx) = 0;

    % right hand side of the mildly nonlinear system
    % inside of the porous medium:
    rhs(1   ,:) = theta(1   ,:) + dt/dz*( Kz(2,:)      - 0          );
    rhs(2:Nz,:) = theta(2:Nz,:) + dt/dz*( Kz(3:Nz+1,:) - Kz(2:Nz,:) );
    % inside of the shallow water layer:
    rhs(Nz+1,:) = theta(Nz+1,:) + dt/dz*( Kz(Nz+2,:) - Kz(Nz+1,:) ) ...
                 -dt/dx*( Hx_tilde(2:Nx+1).*ustar(2:Nx+1) - Hx_tilde(1:Nx).*ustar(1:Nx) );

end




















