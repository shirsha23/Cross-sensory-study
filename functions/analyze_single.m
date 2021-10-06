function [AnaData]=analyze_single_cross_sensory(filename,filename_full_out,data,plotfig)

global is_new_pfit %shir - for using the new pfit
is_new_pfit = 1;

PLOTALL=1;

if nargin <4
    plotfig=1;
end

AnaData=getRightwardChoices_staircase(data);

if plotfig
    filename_fig=filename_full_out(1:end-4);
    if PLOTALL || size(dir(strcat(filename_fig,'.tif')),1)==0 %if the file doesn't already exist
        figure
        set(gcf,'Units','normalized','NumberTitle','off'); 
        set(gcf,'color','w','PaperPositionMode','auto');
        vesL=[0 0 1]; vesR=[0 0 0.3];
        visL=[1 0 0]; visR=[0.3 0 0];
        x = [-16 -8 -4 -2 -1 1 2 4 8 16];
        linewidth=3;
        fontsize=14;
        
        % plots
        try % vestibular
            hold on
            plot(AnaData.xiVespriorL, AnaData.pfitcurveVespriorL,'-','color',vesL,'linewidth',linewidth);
            plot(AnaData.xiVespriorR, AnaData.pfitcurveVespriorR,'-','color',vesR,'linewidth',linewidth);
            scatter(AnaData.dirVespriorL, AnaData.rightVespriorL,AnaData.RepNumVespriorL*5,vesL,'linewidth',linewidth);
            scatter(AnaData.dirVespriorR, AnaData.rightVespriorR,AnaData.RepNumVespriorR*5,vesR,'linewidth',linewidth);
            title('Vestibular test trials','fontsize',fontsize);
            set(gca,'XTick',x,'fontsize',fontsize);
            xlabel('Heading Angle (deg)','fontsize',fontsize);
            y=0:0.5:1;
            ylim([0 1]);
            xlim([-16 16]);
            set(gca,'YTick',y,'fontsize',fontsize);
            ylabel('Rightward Choice Ratio','fontsize',fontsize);
            plot([0 0],[0 1],':','color',[.5 .5 .5],'linewidth',1);
            plot([min(xlim) max(xlim)],[.5 .5],':','color',[.5 .5 .5],'linewidth',1);
            legend('Left Priors','Right Priors','location','southeast');
        end
        
        try % visual
            hold on
            plot(AnaData.xiVisualpriorL, AnaData.pfitcurveVisualpriorL,'-','color',visL,'linewidth',linewidth);
            plot(AnaData.xiVisualpriorR, AnaData.pfitcurveVisualpriorR,'-','color',visR,'linewidth',linewidth);
            scatter(AnaData.dirVisualpriorL, AnaData.rightVisualpriorL,AnaData.RepNumVisualpriorL*5,visL,'linewidth',linewidth);
            scatter(AnaData.dirVisualpriorR, AnaData.rightVisualpriorR,AnaData.RepNumVisualpriorR*5,visR,'linewidth',linewidth);
            title('Visual test trials','fontsize',fontsize);
            set(gca,'XTick',x,'fontsize',fontsize);
            xlabel('Heading Angle (deg)','fontsize',fontsize);
            y=0:0.5:1;
            ylim([0 1]);
            xlim([-16 16]);
            set(gca,'YTick',y,'fontsize',fontsize);
            ylabel('Rightward Choice Ratio','fontsize',fontsize);
            plot([0 0],[0 1],':','color',[.5 .5 .5],'linewidth',1);
            plot([min(xlim) max(xlim)],[.5 .5],':','color',[.5 .5 .5],'linewidth',1);
            legend('Left Priors','Right Priors','location','southeast');
        end
        
        print(gcf,'-dtiff','-r300', filename_fig);
        print(gcf,'-dpdf','-r300', filename_fig);
        savefig(gcf,filename_fig);
        
        if PLOTALL
            close(gcf)
        end
    end
end