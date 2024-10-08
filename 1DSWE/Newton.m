% find h* of the star-region
% h* is the root of G(h*) = Phi_L + Phi_R + uR - uL= 0
function hstar = Newton(hL,hR,uL,uR)

    tol = 1e-10;
    hstar = (hL + hR)/2;
    MaxIter = 100;
    for i=1:MaxIter
        residual = G(hstar,hL,hR,uL,uR);
        if abs(residual) < tol
            break
        end
        dhstar =-residual/dG(hstar,hL,hR,uL,uR);
        d = 1;
        % optimization of the Newton method
        for inner=1:MaxIter
            if (abs(G(hstar+d*dhstar,hL,hR,uL,uR)) >= abs(residual))
                % residual is growing - NOT good
                d = 0.5*d;
            else
                % residual is decreasing - good
                hstar = hstar + d*dhstar;
                break
            end
        end
    end
end