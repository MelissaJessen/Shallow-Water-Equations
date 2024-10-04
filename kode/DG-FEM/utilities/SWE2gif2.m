function SWE2gif(x,ht,hut,b,tt,plotpath,defplotname,FinalTime,N,flux,SL,BC,bottom_geometry,n_frames,save_toggle,hrange,hurange,hubrange)
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
    


for i = 1:length(tt)
    
printcount = printcount + 1;
%     if (printcount == 10*FinalTime||i==length(tt))
    if (printcount == floor(length(tt)/n_frames)||i==length(tt))
        frameidx = frameidx + 1;

        subplot(311)
        plot(x,squeeze(ht(i,:,:))+b)
        hold on
        plot(x,b,'k--')
        hold off
        ylabel('$h+b$','interpreter','latex')
        xlim([xmin,xmax])
        ylim([0,hbmax].*[0.9,1.1])
        grid minor

        subplot(312)
        plot(x,squeeze(ht(i,:,:)))
        ylabel('$h$','interpreter','latex')
        xlim([xmin,xmax])
        ylim([hmin,hmax].*[0.9,1.1])
        grid minor

        subplot(313)
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