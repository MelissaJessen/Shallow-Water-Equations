function [f,fd]=geofun_PREX(he,hk,ak,gravit)

if he<=hk %wave is a rarefaction wave
    a=sqrt(gravit*he);
    f=2./3.*(he*a-hk*ak);
    fd=a;
else %wave is a shock
    ges=sqrt(0.5*gravit*(he+hk));
    f=(he-hk)*ges;
    fd= sqrt(gravit/8)*(3.0*he+hk)/sqrt(he+hk);
end

   