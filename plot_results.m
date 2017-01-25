%plot_results - Plot results for algorithm output
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 18-Mar-2016

%Save GUI Figure to temp_fig
if handles.gui_front==0
    fig_name=['temp_fig/',datestr(datetime),'_aaa_parameters.pdf'];
    warning off
    %Create figure to display parameters etc
    f = figure('Units','characters',...
        'Position',[0 0 150 60],...
        'HandleVisibility','callback',...
        'IntegerHandle','off',...
        'Renderer','painters');
    
    mainPanel = uipanel('Units','characters',...
        'Position',[0 10 90 50],...
        'parent',f);
    botPanel = uipanel('BorderType','etchedin',...
        'Units','characters',...
        'Position',[0 0 150 10],...
        'Parent',f);
    rightPanel = uipanel('BorderType','etchedin',...
        'Units','characters',...
        'Position',[90 10 60 50],...
        'Parent',f);
    
    %Copy parameters
    axesobject_par1=copyobj(handles.values.Children,mainPanel);
    axesobject_fig1=copyobj(handles.motif_image,botPanel);
    set(axesobject_fig1,'Position',[0 0 1 1]);
    panax1 = axes('Units','normal', 'Position', [0 0 1 1], 'Parent', rightPanel);
    set(gca,'Visible','off')
    text(0.2,0.8,'MI-Sim Report','fontsize',24)
    text(0.2,0.7,datestr(datetime),'fontsize',20)
    text(0,0.6,'Equations','fontsize',24)
    text(0,0.4,handles.eqtx,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
    text(0,0.2,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
    
    opts.Width = 20;
    opts.Height = 20;
    opts.format = 'pdf';
    hgexport(f, fig_name, opts);
    close(f)
    handles.gui_front=1;
end
yout=handles.yout;
tout=handles.tout;
fixed_numerical=handles.fn;
%Check that data is available
if exist('yout')==0
    return
end

%For single-point (dynamic) analysis, check which variables to plot
s1_check=get(handles.s1_check,'Value');
x1_check=get(handles.x1_check,'Value');
s2_check=get(handles.s2_check,'Value');
x2_check=get(handles.x2_check,'Value');
s3_check=get(handles.s3_check,'Value');
x3_check=get(handles.x3_check,'Value');
s4_check=get(handles.s4_check,'Value');
s5_check=get(handles.s5_check,'Value');
s6_check=get(handles.s6_check,'Value');
labels=handles.var_names;
checkind=[s1_check,x1_check,s2_check,x2_check,s3_check,x3_check,s4_check,s5_check,s6_check];
fullcheck=sum(checkind);
indices_check = find(checkind);
checkenb=[];
if strcmp(get(handles.s1_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x1_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s2_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x2_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s3_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.x3_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s4_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s5_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
if strcmp(get(handles.s6_check,'Enable'),'on')
    checkenb=[checkenb,1];
else
    checkenb=[checkenb,0];
end
ly=length(yout);
[lfn,sfn]=size(fixed_numerical);
kk=1;
for k=1:length(checkind)
    if checkenb(k)==1
        yout_n(:,k)=yout(:,kk);
        fixed_numerical_n(:,k)=fixed_numerical(:,kk);
        kk=kk+1;
    else
        yout_n(:,k)=zeros(ly,1);
        fixed_numerical_n(:,k)=zeros(lfn,1);
    end
end

%% Solutions plot
%Plots according to routine
switch which_plot
    case 'sol' %Single-point solutions
        switch handles.plotsolution
            
            case 'time' %Time series plot
                try
                    axes(handles.solutionplot); 
                catch
                    subplot(1,1,1)
                    handles.solutionplot=gca;
                end
                %reset the plot
                cla reset
                
                %Variables to plot
                legd=[];
                if s1_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,1),'linewidth',2,'color','b')
                    hold off
                else
                end
                if x1_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,2),'linewidth',2,'color','r')
                    hold off
                else
                end
                if s2_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,3),'linewidth',2,'color',[0 0.5 0])
                    hold off
                else
                end
                if x2_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,4),'linewidth',2,'color',[1 0.6 0])
                    hold off
                else
                end
                if s3_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,5),'linewidth',2,'color',[0.5 0 0.5])
                    hold off
                else
                end
                if x3_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,6),'linewidth',2,'color',[0.302 0.745 0.933])
                    hold off
                else
                end
                if s4_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,7),'linewidth',2,'color',[0 1 1])
                    hold off
                else
                end
                if s5_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,8),'linewidth',2,'color',[1 0.6 0.784])
                    hold off
                else
                end
                if s6_check==1
                    hold on
                    plot(handles.solutionplot,tout,yout_n(:,9),'linewidth',2,'color',[0.314 0.314 0.314])
                    hold off                    
                else
                end
                legd=handles.var_names;
                hold on
                %adds axis labels
                xlabel(handles.solutionplot,'Time (days)')
                ylabel(handles.solutionplot,'Concentration')
                ylimit=get(handles.solutionplot,'Ylim');
                set(handles.solutionplot,'Ylim',[0 ylimit(2)]);
                handles.legend1=legend(legd,'location','best');
                hold off
                
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_dynamics_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj([handles.legend1,handles.solutionplot],newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
                
            case 'deltaG'
                warning off
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                
                plot(handles.tout,handles.dG0_out,'linewidth',2,'color','r')
                hold on
                zr=zeros(1,length(handles.tout));
                plot(handles.tout,zr,'k--')
                xlabel('Time (days)')
                ylabel('\Delta G (kJ/mol)')
                
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_deltaG_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj([handles.solutionplot],newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
                
            case 'inhibition'
                warning off
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                
                plot(handles.tout,handles.I20,'linewidth',2,'color',[0.6 0.2 0.4])
                yl=get(gca,'YLim');
                if max(real(handles.I20))<0
                    mi20=min(real(handles.I20));
                    ma20=0;
                else
                    mi20=0;
                    ma20=max(real(handles.I20));
                end
                set(gca,'YLim',[mi20 ma20])
                xlabel('Time (days)')
                ylabel('Inhibition function')
                axis tight
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_inhibtion_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj([handles.solutionplot],newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
                
            case 'phase' %Phase plot
                
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                
                if fullcheck == 0 || fullcheck == 1
                    msgbox('Too few variables selected, please select 2 or 3 variables','Error','error')
                    return
                elseif fullcheck > 3
                    msgbox('Too many variables selected, please select 2 or 3 variables','Error','error')
                    return
                elseif fullcheck == 2
                    plot(handles.solutionplot,yout_n(:,indices_check(1)),yout_n(:,indices_check(2)),'linewidth',2,'color','b')
                    ylimit=get(handles.solutionplot,'Ylim');
                    set(handles.solutionplot,'Ylim',[0 ylimit(2)]);
                    xlabel(handles.solutionplot,labels(indices_check(1)))
                    ylabel(handles.solutionplot,labels(indices_check(2)))
                elseif fullcheck == 3
                    plot3(handles.solutionplot,yout_n(:,indices_check(1)),yout_n(:,indices_check(2)),yout_n(:,indices_check(3)),'linewidth',2,'color','r')
                    grid on
                    xlabel(handles.solutionplot,labels(indices_check(1)))
                    ylabel(handles.solutionplot,labels(indices_check(2)))
                    zlabel(handles.solutionplot,labels(indices_check(3)))
                end
                
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_phase_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj(handles.solutionplot,newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
            case 'subp' %Subplotting each variable separately
                
                if isfield(handles,'solutionplot') && ishandle(handles.solutionplot)
                    axes(handles.solutionplot)
                elseif isfield(handles,'hsub') 
                   hold on; axes(handles.hsub(1))
                end
                nstat=get(handles.normeig,'Value');

                chk=[s1_check,x1_check,s2_check,x2_check,s3_check,x3_check,s4_check,s5_check,s6_check];
                vars=handles.var_names;
                vc=vars(find(chk));
                
                N = sum(chk);
                
                sa=ceil(N/2); sb=2;

                list=lines(6);
                if handles.clcyc>6
                    list=lines(handles.clcyc);
                end
                
                if nstat==0
                    handles.clcyc=1;
                    for k=1:N
                        handles.hsub(k)=subplot(sa,sb,k);
                        plot(handles.hsub(k),tout,yout_n(:,k),'linewidth',2,'color',list(handles.clcyc,:))
                        ylim=get(gca,'Ylim');
                        handles.hlegend(k)=legend(vc(k)); hold off; set(gca,'Xlim',[0 tout(end)],'Ylim',[0 ylim(2)]);

                        handles.plotaxes_sub{k}=get(gca,'Children');
                            
                    end
                elseif nstat==1
                    handles.clcyc=handles.clcyc+1; hold all
                    for k=1:N
                        
                        axes(handles.hsub(k))
                        hold on
                        
                        plot(handles.hsub(k),tout,yout_n(:,k),'linewidth',2,'color',list(handles.clcyc,:))
                        ylim=get(gca,'Ylim');
                        set(gca,'Xlim',[0 tout(end)],'Ylim',[0 ylim(2)]);

                        handles.plotaxes_sub{k}=get(gca,'Children');
                    end
                end

                %Adds axis labels in centre of subplots
                AxesH    = findobj(handles.hsub, 'Type', 'Axes');
                YLabelHC = get(AxesH, 'YLabel'); XLabelHC = get(AxesH,'XLabel');
                YLabelH  = [YLabelHC{:}]; XLabelH  = [XLabelHC{:}];
                set(YLabelH, 'String', 'Conc.')
                set(XLabelH, 'String', 'Time (d)')
                
                %suplabel('Time (days)');
                %suplabel('Concentration (kgCOD m^{-3})','y');
                
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_subplot_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj(handles.solutionplot,newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
            case 'eig' %Eigenvalue plot
                axes(handles.solutionplot);
                %reset the plot
                cla reset
                nstat=get(handles.normeig,'Value');
                eig=handles.eigen';
                %Get real and imaginary parts
                if nstat==0
                    re = real(eig);
                    ime = imag(eig);
                    xlabel('Real')
                    ylabel('Img')
                else
                    re=normr(real(eig));
                    ime=imag(eig);
                end
                
                [m,n]=size(eig);
                col=num2cell(jet(n),2);
                markers = {'o','s','d','v','x','+','*','^','.','>','<','p','h'};
                mrk=markers(1:n)';
                lh=plot(re,ime,'.','markersize',10);
                set(lh,{'marker'},mrk,{'markerfacecolor'},col,'markeredgecolor',[0,0,0]);

                if nstat==0
                    xlabel('Real')
                    ylabel('Img')
                else
                    xlabel('Norm(Real)')
                    ylabel('Img')
                end
                
               %Legend
               for k=1:n
                   lb(k)=cellstr(['FP',num2str(k)]);
               end
               legend(lb,'location','Best')
               xl=get(gca,'Xlim');
               yl=get(gca,'Ylim');
               hold on
               hx=line(xl,[0 0]);
               hy=line([0 0],yl);
               set(hx,'linestyle','--','color','k')
               set(hy,'linestyle','--','color','k')
               
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_eigenvalue_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj(handles.solutionplot,newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
        end

    case 'traj'  %Trajectory plot
        pben=get(handles.plot_button,'enable');
        if strcmp(pben,'on')==0
            return
        end
        axes(handles.trajectoryplot);
        
        switch handles.plotdim
            case 'two' %2D plot
                
                %check if plot is to be overlayed
                multi=get(handles.overlay,'Value');
                if multi==0
                    cla reset
                else
                end
                
                %Check correct number of variables are selected
                if fullcheck ~=2
                    msgbox('Too few variables selected, please select 2 variables','Error','error')
                    return
                else
                    X_value=yout_n(:,indices_check(1));
                    Y_value=yout_n(:,indices_check(2));
                    hold on
                    %plot the fixed points
                    z=length(fixed_numerical(:,1));
                    list=hsv(z);
                    for i=1:z
                        plot(handles.trajectoryplot,fixed_numerical_n(i,indices_check(1)),fixed_numerical_n(i,indices_check(2)),'s','color',list(i,:),'MarkerFaceColor',list(i,:))
                    end
                    %plot the trajectory
                    plot(handles.trajectoryplot,X_value,Y_value,'linewidth',2)
                    xlabel(handles.trajectoryplot,labels(indices_check(1),:))
                    ylabel(handles.trajectoryplot,labels(indices_check(2),:))
                    %plot the start point (green dot)
                    plot(handles.trajectoryplot,yout_n(1,indices_check(1)),yout_n(1,indices_check(2)),'o','color','g','MarkerFaceColor','g')
                    %plot the end point (red dot)
                    plot(handles.trajectoryplot,yout_n(end,indices_check(1)),yout_n(end,indices_check(2)),'o','color','r','MarkerSize',10,'linewidth',2)
                    hold off
                end
                set(handles.overlay,'enable','on')
                
                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_2d_trajectory_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
            case 'three' %3D plot
                %if user wants to overlay multiple trajectories onto the plot of
                %trajectories
                multi=get(handles.overlay,'Value');
                if multi==0
                    cla reset
                else
                    hold on
                end
                
                %Check correct number of variables are selected
                if fullcheck ~=3
                    msgbox('Too few variables selected, please select 3 variables','Error','error')
                    return
                else
                    X_value=yout_n(:,indices_check(1));
                    Y_value=yout_n(:,indices_check(2));
                    Z_value=yout_n(:,indices_check(3));
                    
                    %plot the fixed points
                     z=length(fixed_numerical(:,1));
                    list=hsv(z);
                    for i=1:z
                        plot3(handles.trajectoryplot,squeeze(fixed_numerical_n(i,indices_check(1))),squeeze(fixed_numerical_n(i,indices_check(2))),squeeze(fixed_numerical_n(i,indices_check(3))),'s','color',list(i,:),'MarkerFaceColor',list(i,:))
                        hold on
                    end
                    %plot the trajectory
                    plot3(handles.trajectoryplot,X_value,Y_value,Z_value,'linewidth',2)
                    xlabel(handles.trajectoryplot,labels(indices_check(1)))
                    ylabel(handles.trajectoryplot,labels(indices_check(2)))
                    zlabel(handles.trajectoryplot,labels(indices_check(3)))
                    %plot the start point (green dot)
                    plot3(handles.trajectoryplot,yout_n(1,indices_check(1)),yout_n(1,indices_check(2)),yout_n(1,indices_check(3)),'o','color','g','MarkerFaceColor','g')
                    %plot the end point (red dot)
                    plot3(handles.trajectoryplot,yout_n(end,indices_check(1)),yout_n(end,indices_check(2)),yout_n(end,indices_check(3)),'o','color','r','MarkerSize',10,'linewidth',2)
                    grid on
                    hold off
                end
                set(handles.overlay,'enable','on')

                %Save Figure to temp_fig
                fig_name=['temp_fig/',datestr(datetime),'_3d_trajectory_plot.pdf'];
                newfig_a=figure;
                axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
                hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
                close(newfig_a)
        end
        
end
handles.xlabelhandle=get(gca,'Xlabel');
handles.ylabelhandle=get(gca,'Ylabel');
handles.zlabelhandle=get(gca,'Zlabel');
handles.plothandle=get(gca,'Children');
guidata(hObject,handles);


