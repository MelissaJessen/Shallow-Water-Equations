function [SL,SR]=hll(hL,hR,uL,uR,aL,aR,gravit)
    

hstar=(((aL+aR)/2+(uL-uR)/4)^2)/gravit;

% hS2 = (0.75/sqrt(gravit)*(uL*hL-uR*hR) + 0.5*(hL^1.5+hR^1.5))^2;
% hstar = hS2^(1/3);
    
    if hstar > hL
        qL=sqrt(0.5*(hstar+hL)*hstar/(hL^2));
    else
        qL=1;
    end
    if hstar > hR
        qR=sqrt(0.5*(hstar+hR)*hstar/(hR^2));
    else
        qR=1;
    end
    SL=uL   - aL*qL;
    SR=uR + aR*qR;
end