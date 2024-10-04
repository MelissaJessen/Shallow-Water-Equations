     %%% Initializes variables for Advection DG-FEM solver

[Nv, VX, K, EToV] = MeshGen1D(xmin,xmax,K);
Np        = N+1;
Nfp       = 1;
Nfaces    = 2;
r         = JacobiGL(0,0,N);
V         = Vandermonde1D(N, r); 
invV      = inv(V);
Dr        = Dmatrix1D(N, r, V);
Emat      = zeros(Np,Nfaces*Nfp);
Emat(1,1) = 1.0; Emat(Np,2) = 1.0;
Minv      = V*V';
M         = inv(Minv);

                         % a per element basis otherwise

% surfint   = V*(V'*Emat); % Take-out elements associated with boudnaries
surfint   = Minv*Emat;

va        = EToV(:,1)'; vb = EToV(:,2)';
x         = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));
[rx,J]    = GeometricFactors1D(x,Dr);
Mk        = 0.5*J(1).*M; % Only true for uniform elements, Mk needs to be computed on eacg element

fmask1    = find( abs(r+1) < NODETOL)'; 
fmask2    = find( abs(r-1) < NODETOL)';
Fmask     = [fmask1;fmask2]';
Fx        = x(Fmask(:), :);
[nx]      = Normals1D(Nfp,Nfaces,K);
Fscale    = 1./(J(Fmask,:));
[EToE, EToF] = Connect1D(EToV);
[vmapM, vmapP, vmapB, mapB,mapI,mapO,vmapI,vmapO] = BuildMaps1D(x,Np,Nfp,Nfaces,K,Fmask,EToE,EToF,NODETOL);
if ~periodic
    param.periodic = periodic;
else
    param.periodic = 1;
    vmapP(1) = vmapM(end);
    vmapP(end) = vmapM(1);
end
