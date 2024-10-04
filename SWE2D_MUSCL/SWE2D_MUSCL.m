% MUSCL scheme for the 2D SWE 
clear;
clc;

global g Nx Ny dx dy;

% physical parameters
g  = 9.81;       % acceleration due to gravity
xL =-1; xR = 1;
yL =-1; yR = 1;
t  = 0;
tend = 1;

% numerical parameters
Nx = 64;
Ny = Nx;
dx = (xR - xL)/Nx;
dy = (yR - yL)/Ny;
x  = linspace(xL+dx/2,xR-dx/2,Nx);
y  = linspace(yL+dy/2,yR-dy/2,Ny);
CFL = 0.9;  % CFL number < 1

% Initial conditions
Q = zeros(3,Nx,Ny);
Qnew = Q;
for i=1:Nx
    for j=1:Ny
        r = sqrt( (x(i) - 0.9 )^2 + (y(j) - 0.9 )^2);
        if ( r < 0.5 )
            Q(1,i,j) = 2;
        else
            Q(1,i,j) = 0.2;
        end
    end
end

surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp');

% allocate auxiliary arrays for the MUSCL part
slopeX = zeros(3,Nx,Ny);
slopeY = zeros(3,Nx,Ny);
Qxm = zeros(3,Nx,Ny);
Qxp = zeros(3,Nx,Ny);
Qym = zeros(3,Nx,Ny);
Qyp = zeros(3,Nx,Ny);

% time loop
for n=1:100000
    sx = Lambdax(Q);    % compute characteristic velocities in the x-directio
    sy = Lambday(Q);    % compute characteristic velocities in the y-directio
    ax = max( max ( max( abs(sx) ) ) );
    ay = max( max ( max( abs(sy) ) ) );
    dt = CFL/(ax/dx + ay/dy);
    if (t + dt > tend)
        dt = tend - t;
    end
    if (t >= tend)
        break
    end

    % MUSCL part:
    % slope limiter
    slopeX(:,2:Nx-1,:) = minmodarray( Q(:,2:Nx-1,:) - Q(:,1:Nx-2,:) , Q(:,3:Nx,:) - Q(:,2:Nx-1,:) );
    slopeY(:,:,2:Ny-1) = minmodarray( Q(:,:,2:Ny-1) - Q(:,:,1:Ny-2) , Q(:,:,3:Ny) - Q(:,:,2:Ny-1) );
    % extrapolation of the data toward the left and right faces of the cell
    Qxm = Q - 0.5*slopeX;
    Qxp = Q + 0.5*slopeX;
    Qym = Q - 0.5*slopeY;
    Qyp = Q + 0.5*slopeY;
    % extrapolation in time
    Q_t = - ( FluxX(Qxp) - FluxX(Qxm) )/dx - ( FluxY(Qyp) - FluxY(Qym) )/dy;
    Qxm = Qxm + 0.5*dt*Q_t;
    Qxp = Qxp + 0.5*dt*Q_t;
    Qym = Qym + 0.5*dt*Q_t;
    Qyp = Qyp + 0.5*dt*Q_t;

    % recompute the physical fluxes and characteristic velocities (this we need for the Rusanov flux)
    fxm = FluxX(Qxm); fxp = FluxX(Qxp);
    gym = FluxY(Qym); gyp = FluxY(Qyp);
    
    sxm = Lambdax(Qxm); sxp = Lambdax(Qxp);
    sym = Lambday(Qym); syp = Lambday(Qyp);
    
    % space loop
    dtdx = dt/dx;
    dtdy = dt/dy;    
    for i = 1:Nx
        for j = 1:Ny
            % numerical fluxes in X
            if     ( i == 1)
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(2) =-Qghost(2);  % flip the sign of the x momentum
                Fm = Rusanov(Qghost      ,Qxm(:,i  ,j),FluxX(Qghost),fxm(:,i  ,j),Lambdax(Qghost),sxm(:,i  ,j));
                Fp = Rusanov(Qxp(:,i  ,j),Qxm(:,i+1,j),fxp(:,i  ,j) ,fxm(:,i+1,j),sxp(:,i  ,j)   ,sxm(:,i+1,j));
            elseif (i == Nx)
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(2) =-Qghost(2);  % flip the sign of the x momentum
                Fm = Rusanov(Qxp(:,i-1,j),Qxm(:,i  ,j),fxp(:,i-1,j),fxm(:,i  ,j) ,sxp(:,i-1,j),sxm(:,i ,j)    );
                Fp = Rusanov(Qxp(:,i  ,j),Qghost    ,fxp(:,i  ,j),FluxX(Qghost),sxp(:,i  ,j),Lambdax(Qghost));
            else
                Fm = Rusanov(Qxp(:,i-1,j),Qxm(:,i  ,j),fxp(:,i-1,j),fxm(:,i  ,j),sxp(:,i-1,j),sxm(:,i ,j));
                Fp = Rusanov(Qxp(:,i  ,j),Qxm(:,i+1,j),fxp(:,i  ,j),fxm(:,i+1,j),sxp(:,i  ,j),sxm(:,i+1,j));
            end
            % numerical fluxes in Y
            if     ( j == 1  )
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(3) =-Qghost(3);  % flip the sign of the x momentum
                Gm = Rusanov(Qghost      ,Qym(:,i,j  ),FluxY(Qghost),gym(:,i,j  ),Lambday(Qghost),sym(:,i,j  ));
                Gp = Rusanov(Qyp(:,i,j  ),Qym(:,i,j+1),gyp(:,i,j  ) ,gym(:,i,j+1),syp(:,i,j  )   ,sym(:,i,j+1));
            elseif ( j == Ny )
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(3) =-Qghost(3);  % flip the sign of the x momentum
                Gm = Rusanov(Qyp(:,i,j-1),Qym(:,i,j  ),gyp(:,i,j-1),gym(:,i,j  ) ,syp(:,i,j-1),sym(:,i,j  )   );
                Gp = Rusanov(Qyp(:,i,j  ),Qghost    ,gyp(:,i,j  ),FluxY(Qghost),syp(:,i,j  ),Lambday(Qghost));
            else
                Gm = Rusanov(Qyp(:,i,j-1),Qym(:,i,j  ),gyp(:,i,j-1),gym(:,i,j  ),syp(:,i,j-1),sym(:,i,j  ));
                Gp = Rusanov(Qyp(:,i,j  ),Qym(:,i,j+1),gyp(:,i,j  ),gym(:,i,j+1),syp(:,i,j  ),sym(:,i,j+1));
            end
            % stadard finite volume update
            Qnew(:,i,j) = Q(:,i,j) - dtdx*( Fp - Fm ) - dtdy*( Gp - Gm ) ;
        end
    end
    % update time
    t = t + dt;
    Q = Qnew;

    % plot the solution
    subplot(2,1,1)
    surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp');
    title(strcat('t = ',num2str(t)))
    xlabel('x')
    ylabel('y')
    zlabel('h')
    axis([xL xR yL yR 0 2])
    camlight;

    subplot(2,1,2)
    plot(x,Q(1,:,Ny/2))
    xlabel('x')
    ylabel('h(t,x,0)')

    pause(0.001)     
end







