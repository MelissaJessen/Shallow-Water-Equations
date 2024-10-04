function [nx] = Normals1D(Nfp,Nfaces,K)

% function [nx] = Normals1D
% Purpose : Compute outward pointing normals at elements faces

nx = zeros(Nfp*Nfaces, K); 

% Define outward normals
nx(1, :) = -1.0; nx(2, :) = 1.0;
return
