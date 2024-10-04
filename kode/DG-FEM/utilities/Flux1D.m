function du = Flux1D(um,up,nx,alpha)

[avg,jump] = AvgJump(um,up,nx);
du = avg + 0.5*(1-alpha)*jump;
end