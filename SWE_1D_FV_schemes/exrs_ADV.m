function [hstar,ustar]=exrs_ADV(hL,hR,uL,uR)

R_h = hR/hL;
R_U = uL/uR;

if (uL>=0 && uR>=0 )
    %solution in the star region
    hstar = hR;
    ustar = uL*sqrt(hL/hR);
elseif (uL<0 && uR<0 )
    %solution in the star region
    hstar = hL;
    ustar = -abs(uR)*sqrt(hR/hL);
elseif (uL>0 && uR<0 )
    if (abs(R_U) < sqrt(R_h)) %left wave is a shock wave
        %solution in the star region
        hstar = hL;
        ustar = -abs(uR)*sqrt(hR/hL);
    else %right wave is a shock wave
        %solution in the star region ??Non funziona!!!!!
        hstar = hR;
        ustar = uL*sqrt(hL/hR);
    end
elseif (uL<0 && uR>0 ) % transonic rarefaction
    %solution in the star region
    hstar = hL;% or hR It works well with both choices
    ustar = 0;
end

end

    
       





















