function [d,u] = drybed(d,u,param)
chalen = param.chalen;
mcells = param.mcells;
dr = param.dr;
dl = param.dl;

% fld,fld,fr,frd i dont think these need to be defined as they are computed
% each step in geofun;
gate = param.gate;
timout = param.timout;


for i = 1:mcells
    xcoord = i*chalen/mcells - gate;
    s      = xcoord/timout;
%     di = d(i);
%     ui = u(i);
    if dl <= 0
        % left state is dry
        [dsam,usam] = samlef(s,param);
    else
        if(dr<=0)
            % right state is dry
            [dsam,usam] = samrig(s,param);
        else
            %middle state is dry
            [dsam,usam] = sammid(s,param);
        end
    end
    d(i) = dsam;
    u(i) = usam;
end

end