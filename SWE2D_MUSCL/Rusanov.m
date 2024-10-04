% Rusanov numerical flux
function F = Rusanov(QL,QR,FL,FR,sL,sR)

    smax = max( max( abs(sL) ),max( abs(sR) ) );

    F = 0.5*( FR + FL ) - 0.5*smax*(QR - QL);

end