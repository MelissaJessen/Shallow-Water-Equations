% compute the exact solution of the Riemann problem for the SW equations
% QL - is the left constant state
% QR - is the right constant state
% lambda = x/t
function Q = ExactRP(QL,QR,lambda)

    global g;

% compute the star solution Q* 
hL = QL(1);
hR = QR(1);
uL = QL(2)/hL;
uR = QR(2)/hR;
hstar = Newton(hL,hR,uL,uR);
ustar = uL - Phi(hstar,hL);

% sampling of the left wave
if (lambda <= ustar) % sampling of the left wave
    if (hstar > hL) % then the left wave is a shock
        s = uL - sqrt(0.5*g*hstar/hL*(hstar + hL));
        if (lambda <= s)
            Q = [hL; hL*uL];
        else
            Q = [hstar; hstar*ustar];
        end
    else % then the left wave is a rarefaction
        urw_head = uL    - sqrt(g*hL   ); % velocity of the head of the rarefaction wave
        urw_tail = ustar - sqrt(g*hstar); % velocity of the tail of the rarefaction wave
        if (lambda <= urw_head)
            % x_{i+1/2} is in the QL region
            Q = [hL; hL*uL];    
        elseif(lambda >= urw_tail)
            % x_{i+1/2} is in the Qstar region
            Q = [hstar; hstar*ustar];
        else
            % x_{i+1/2} is in the rarefaction wave
            h = ((uL + 2*sqrt(g*hL) - lambda)/3)^2/g;
            Q = [h; h*(lambda + sqrt(g*h))];
        end
    end
%% *** sampling of the right wave ***
else        
    if (hstar > hR) % then the right wave is a shock
        s = uR + sqrt(0.5*g*hstar/hR*(hstar + hR));
        if (lambda <= s)
            Q = [hstar; hstar*ustar];
        else
            Q = [hR; hR*uR];
        end
    else % then the left wave is a rarefaction
        urw_head = uR    + sqrt(g*hR   ); % velocioty of the head of the rarefaction wave
        urw_tail = ustar + sqrt(g*hstar); % velocity of the tail of the rarefaction wave
        if (lambda >= urw_head)
            % x_{i+1/2} is in the QL region
            Q = [hR; hR*uR];    
        elseif(lambda <= urw_tail)
            % x_{i+1/2} is in the Qstar region
            Q = [hstar; hstar*ustar];
        else
            % x_{i+1/2} is in the rarefaction wave
            h = ((-uR + 2*sqrt(g*hR) + lambda)/3)^2/g;
            Q = [h; h*(lambda - sqrt(g*h))];
        end
    end
end

end