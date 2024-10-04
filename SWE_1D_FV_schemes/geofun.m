function [f,fd]=geofun(he,hk,ck,gravit)

if he<=hk %wave is a rarefaction wave
    c=sqrt(gravit*he);
    f=2*(c-ck);
    fd=gravit/c;
else %wave is a shock
    ges=sqrt(0.5*gravit*(he+hk)/(he*hk));
    f=(he-hk)*ges;
    fd=ges-0.25*gravit*(he-hk)/(ges*he^2);
end

    