% MUSCL scheme for the 2D SWE 
clear;
clc;

global g Nx Ny dx dy;

% physical parameters
g  = 9.81;       % acceleration due to gravity
xL =-1; xR = 1;
yL =-1; yR = 1;
t  = 0;
tend = 0.1;

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
        r = sqrt( (x(i) - 0 )^2 + (y(j) - 0 )^2);
        if ( r < 0.5 )
            Q(1,i,j) = 2;
        else
            Q(1,i,j) = 0.2;
        end
    end
end

surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp');

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

    fxm = FluxX(Q); fxp = FluxX(Q);
    gxm = FluxY(Q); gxp = FluxY(Q);
    
    sxm = Lambdax(Q); sxp = Lambdax(Q);
    sym = Lambday(Q); syp = Lambday(Q);
    
    % space loop
    dtdx = dt/dx;
    dtdy = dt/dy;    
    for i = 1:Nx
        for j = 1:Ny
            % numerical fluxes in X
            if     ( i == 1)
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(2) =-Qghost(2);  % flip the sign of the x momentum
                Fm = Rusanov(Qghost    ,Q(:,i  ,j),FluxX(Qghost),fxm(:,i  ,j),Lambdax(Qghost),sxm(:,i  ,j));
                Fp = Rusanov(Q(:,i  ,j),Q(:,i+1,j),fxp(:,i  ,j) ,fxm(:,i+1,j),sxp(:,i  ,j)   ,sxm(:,i+1,j));
            elseif (i == Nx)
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(2) =-Qghost(2);  % flip the sign of the x momentum
                Fm = Rusanov(Q(:,i-1,j),Q(:,i  ,j),fxp(:,i-1,j),fxm(:,i  ,j) ,sxp(:,i-1,j),sxm(:,i ,j)    );
                Fp = Rusanov(Q(:,i  ,j),Qghost    ,fxp(:,i  ,j),FluxX(Qghost),sxp(:,i  ,j),Lambdax(Qghost));
            else
                Fm = Rusanov(Q(:,i-1,j),Q(:,i  ,j),fxp(:,i-1,j),fxm(:,i  ,j),sxp(:,i-1,j),sxm(:,i ,j));
                Fp = Rusanov(Q(:,i  ,j),Q(:,i+1,j),fxp(:,i  ,j),fxm(:,i+1,j),sxp(:,i  ,j),sxm(:,i+1,j));
            end
            % numerical fluxes in Y
            if     ( j == 1  )
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(3) =-Qghost(3);  % flip the sign of the x momentum
                Gm = Rusanov(Qghost    ,Q(:,i,j  ),FluxY(Qghost),gxm(:,i,j  ),Lambday(Qghost),sym(:,i,j  ));
                Gp = Rusanov(Q(:,i,j  ),Q(:,i,j+1),gxp(:,i,j  ) ,gxm(:,i,j+1),syp(:,i,j  )   ,sym(:,i,j+1));
            elseif ( j == Ny )
                Qghost    = Q(:,i,j);   % create a ghost state 
                Qghost(3) =-Qghost(3);  % flip the sign of the x momentum
                Gm = Rusanov(Q(:,i,j-1),Q(:,i,j  ),gxp(:,i,j-1),gxm(:,i,j  ) ,syp(:,i,j-1),sym(:,i,j  )   );
                Gp = Rusanov(Q(:,i,j  ),Qghost    ,gxp(:,i,j  ),FluxY(Qghost),syp(:,i,j  ),Lambday(Qghost));
            else
                Gm = Rusanov(Q(:,i,j-1),Q(:,i,j  ),gxp(:,i,j-1),gxm(:,i,j  ),syp(:,i,j-1),sym(:,i,j  ));
                Gp = Rusanov(Q(:,i,j  ),Q(:,i,j+1),gxp(:,i,j  ),gxm(:,i,j+1),syp(:,i,j  ),sym(:,i,j+1));
            end
            % stadard finite volume update
            Qnew(:,i,j) = Q(:,i,j) - dtdx*( Fp - Fm ) - dtdy*( Gp - Gm ) ;
        end
    end
    % update time
    t = t + dt;
    Q = Qnew;

    % plot the solution
    surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp');
    title(strcat('t = ',num2str(t)))
    xlabel('x')
    ylabel('y')
    zlabel('h')
    axis([xL xR yL yR 0 2])
    
    pause(0.001)     
end

plot(x,Q(1,:,Ny/2))






