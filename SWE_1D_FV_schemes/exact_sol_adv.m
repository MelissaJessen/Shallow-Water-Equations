function [hexact,uexact,coexact,x] = exact_sol_adv(mcells,hL,hR,uL,uR,coL,coR,gravit,time,xmin,xmax)
                                  

%mcells  =5000;
tol     =1.E-12;

h  = zeros(1,mcells+1);
u  = zeros(1,mcells+1);
co = zeros(1,mcells+1);

%hexact=zeros(1,mcells+1);
%uexact=zeros(1,mcells+1);

% compute ratios
R_h = hR/hL;
if (abs(uR)>tol)
    R_U = uL/uR;
else
    R_U = 0;
end


%Hcrit=(uR-uL)-2*(cL+cR); %depth positivity condition
%scrsz=get(0,'ScreenSize');
%f1=figure('Position',[scrsz(3)/4 10 3*scrsz(4)/4 scrsz(3)]);
if (uL>=0 && uR>=0 )
    %solution in the star region
    hS = hR;
    uS = uL*sqrt(hL/hR);
    if R_U > sqrt(R_h) %right wave is a shock wave
        sR = uR + uS; % shock speed
        right = 1 ;
    else % right wave is a rarefaction
        shR = 2*uR;
        stR = 2*uS;
        right = 0;
    end
elseif (uL<0 && uR<0 )
    %solution in the star region
    hS = hL;
    uS = -abs(uR)*sqrt(hR/hL);
    if R_U < sqrt(R_h) %left wave is a shock wave
        sL = uL + uS; % shock speed
        left = 1 ;
    else % left wave is a rarefaction
        shL = 2*uL;
        stL = 2*uS;
        left = 0;
    end
elseif (uL>0 && uR<0 )
    if (abs(R_U) < sqrt(R_h)) %left wave is a shock wave
        %solution in the star region
        hS = hL;
        uS = -abs(uR)*sqrt(hR/hL);
        sL = uL + uS; % shock speed
        left = 1 ;
    else %right wave is a shock wave
        %solution in the star region
        hS = hR;
        uS = uL*sqrt(hL/hR);
        sR = uR + uS; % shock speed
        right = 1 ;
    end
elseif (uL<0 && uR>0 )
    %solution in the star region
    %hS = ??;
    uS = 0;
% transonic rarefaction
        shL = 2*uL;
        stR = 2*uR;
%        left = 0;   
end
   
    
    dx=(xmax-xmin)/mcells;
    x=xmin:dx:xmax;
    
    S = x/time;
    for i=1:mcells+1
        if uS< -tol
            %co(i) = coL;
            if left==1
                if S(i)<sL
                    h(i)=hL;
                    u(i)=uL;
                elseif (S(i) <=0)
                    h(i)=hS;
                    u(i)=uS;
                else
                    h(i)=hR;
                    u(i)=uR;
                end
            else
                if S(i)<shL
                    h(i)=hL;
                    u(i)=uL;
                else
                    if S(i)<stL
                        u(i)=S(i)/2;
                        h(i)=hL;
                    elseif (S(i) <=0)
                        h(i)=hS;
                        u(i)=uS;
                    else
                        h(i)=hR;
                        u(i)=uR;
                    end
                end
            end
            
        elseif (uS > tol)
            %                co(i) = coR;
            if right==1 %right shock
                if S(i)>sR
                    h(i)=hR;
                    u(i)=uR;
                elseif (S(i) >= 0)
                    h(i)=hS;
                    u(i)=uS;
                else
                    h(i)=hL;
                    u(i)=uL;
                end
            else %right rarefaction
                if S(i)>shR
                    h(i)=hR;
                    u(i)=uR;
                else
                    if S(i)>stR
                        u(i)=S(i)/2;
                        h(i) = hR;
                    elseif (S(i) >= 0)
                        h(i)=hS;
                        u(i)=uS;
                    else
                        h(i)=hL;
                        u(i)=uL;
                    end
                end
            end
        else %(uS=0: sonic rarefaction)
            if S(i)<shL
                h(i)=hL;
                u(i)=uL;
            elseif S(i) >stR
                h(i)=hR;
                u(i)=uR;
            else
                u(i) = S(i)/2;
                if (S(i) < 0.0)
                    h(i) = hL;
                else
                    h(i) = hR;
                end
            end
             
            
        end
        %         end
        
    end
    
    
for i=1:mcells+1
    hexact(i)=h(i);
    uexact(i)=u(i);
    coexact(i)=co(i);
    %x(i)=x(i)-25;
end

end

