function SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,SL,BC,bottom_geometry,n_frames,hrange,hrange,save_toggle)
% close all
% clear mov;
figure
printcount = 0;
if save_toggle
filename  = PlotName(plotpath,defplotname,FinalTime,N,flux,SL,BC,bottom_geometry,'.gif');
end
frameidx = 0;
% printfactor = 100;
xmin = min(min(x));
xmax = max(max(x));

% hmin = min(ht(:));
% hmax = max(ht(:));
% hbmin = min(b(:)) + hmin;
% hbmax = max(b(:)) + hmax;
% humin = min(hut(:));
% humax = max(hut(:));




if length(hrange) == 0
    hmin = min(real(ht(:)));
    hmax = max(real(ht(:)));
    humin = min(real(b(:))) + hmin;
    humax = max(real(b(:))) + hmax;    
else
    hmin  = hrange(1);
    hmax  = hrange(2);
    humin = hurange(1);
    humax = hurange(2);
end
    
    


% if ~(all([isreal(hbmin),isreal(hbmax),isreal(hmin),isreal(hmax),isreal(hubmin),isreal(humax)]))
% hmin = real(hmin);hmax = abs(hmax);hbmin = abs(hbmin);hbmax = abs(hbmax);humin = abs(humin);humax = abs(humax);

% if all(b(:)==0)
for i = 1:length(tt)

    printcount = printcount + 1;
%     if (printcount == 10*FinalTime||i==length(tt))
    if (printcount == floor(length(tt)/n_frames)||i==length(tt))
        frameidx = frameidx + 1;

        subplot(211)
        plot(x,squeeze(ht(i,:,:))+b)
        hold on
        plot(x,b,'k--')
        hold off
        ylabel('$h$','interpreter','latex')
        xlim([xmin,xmax])
        ylim([0,hbmax].*[0.9,1.1])
        grid minor

        subplot(212)
        plot(x,squeeze(hut(i,:,:)))
        ylabel('$hu$','interpreter','latex')
        grid minor
        xlim([xmin,xmax])
        ylim([humin,humax].*[0.9,1.1])
        sgtitle(['T = ',num2str(tt(i)),'$s$'],'interpreter','latex');
        set(gcf,'color','w')

        drawnow
        mov(frameidx) = getframe(gcf);
        printcount = 0;
        pause(0.01)
    end
end

    
if save_toggle
movie2gif(mov, filename, 'LoopCount', 3, 'DelayTime', 0)
end
end