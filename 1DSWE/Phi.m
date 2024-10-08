function F = Phi(hstar,hi)
global g;

    if (hstar <= hi)
        % shock wave
        F = 2*( sqrt(g*hstar) - sqrt(g*hi) );
    else
        % rarefaction wave
        F = sqrt(0.5*g*(hstar + hi)/(hstar*hi))*(hstar - hi);
    end
end