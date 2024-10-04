function [d,u] = samwet(s,ds,us,dl,ul,dr,ur,cl,cr,cs,param)
g = param.g;

if s<=us
    % sample left wave
    if ds >= dl
        %%%%
        % left shock
        %%%%
        ql = sqrt((ds+dl)*ds/(2.0*dl*dl));
        sl = ul - cl*ql;
        if s <=sl
            % sample point lies to the left of the shock
            d = dl;
            u = ul;
        else
            % sample point lies to the right of the shock
            d = ds;
            u = us;
        end
    else
        %%%
        % left rarefaction
        %%%
        shl = ul-cl;
        if s<=shl
            % sample point lies to the right of the rarefaction
            d = dl;
            u = ul;
        else
            stl = us - cs;
            if s <= stl
                % sample point lies inside the rarefaction
                u = (ul + 2 * cl + 2 * s)/3.0;
                c = (ul + 2 * cl -s)/3.0;
                d = c^2/g;
            else
                % sample point lies in the STAR region
                d = ds;
                u = us;
            end
        end
    end
else
    %%%
    % sample right wave
    %%%
     if ds >= dr
         % right shock
         qr = sqrt((ds+dr)*ds/(2*dr^2));
         sr = ur + cr*qr;
%          'test1'
         if s>=sr
             % sample point lies to the right of the shock
             d = dr;
             u = ur;
%              'test2'
         else
             % sample point lies to the left of the shock
             d = ds;
             u = us;
%              'test3'
         end
     else
         % right rarefaction
         shr = ur + cr;
         if s >= shr
             % sample point lies to the right of the rarefaction
             d =  dr;
             u = ur;
%              'test4'
         else
             str = us + cs;
             if s>=str
                % sample ppoint lies inside the rarefaction
                u = (ur-2*cr + 2*s)/3;
%                 'test5'
                c = (-ur + 2*cr + s)/3;
                d = c^2/g;
             else
                 % sample point lies inside the star region
%                 'test6'
                 d = ds;
                 u = us;
             end
         end
     end
end

                
          