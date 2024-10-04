% minmod function required for the limiting procedure 

function c = minmodarray(a,b)

    [nVar,Nx,Ny] = size(a);

    c = zeros(nVar,Nx,Ny);

    for i = 1:nVar
        logic = (a(i,:,:).*b(i,:,:) > 0 & abs(a(i,:,:)) < abs(b(i,:,:)));
        c(i,:,:) = a(i,:,:).*logic;
        logic = (a(i,:,:).*b(i,:,:) > 0 & abs(a(i,:,:)) > abs(b(i,:,:)));
        c(i,:,:) = c(i,:,:) + b(i,:,:).*logic;
    end

end