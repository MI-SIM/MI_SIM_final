function motif_out=model_sel(handles)
%model_sel - initialisation and parameterisation of motifs
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% February 2015; Last revision: 1-Apr-2016

set(handles.infotable,'Visible','off');
set(handles.stabtable,'Visible','off');
switch handles.motif_name
    case 'Commensalism'
        set(handles.model_label,'String','4 ODE: Food Chain');
        set(handles.commensalism,'Checked','on'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='fc';
        axes(handles.motif_image)
        imshow('motifs_images/fc.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','off'); set(handles.S3_0,'enable','off')
        set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','on')
        set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
        set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')

        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','off','String',0);
        set(handles.s3in_in,'Enable','off','String',0);
        set(handles.ki2_in,'String',0);
        set(handles.ks1_in,'String',0.3);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','on','String',35);
        set(handles.y2_in,'Enable','on','String',0.06);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',2.5e-5);
        set(handles.km3_in,'Enable','off','String',0);
        set(handles.y3_in,'Enable','off','String',0);
        set(handles.kdec3_in,'Enable','off','String',0);
        set(handles.ks3_in,'Enable','off','String',0);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','on','String',0.43);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
        set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
        set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
        set(handles.gam_text1,'Enable','off');set(handles.gam_text2,'Enable','off')
        
        handles.plotsolution='time'; handles.plotdim='two';
        
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma_0(1-Y_1)f_1X_1 - f_2X_2$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
        
        %Functions depend on growth model
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        
        handles.params_sim=strvcat('S1in','D','km1','KS1','Y1','km2','KS2','Y2','kdec1','kdec2');
        handles.var_names=strvcat('S1','X1','S2','X2');
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
        cla(handles.txax)
        axes(handles.txax)
        
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        text(0.7,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'Competition'
        set(handles.model_label,'String','3 ODE: Substrate Competition');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','on');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='sc';
        axes(handles.motif_image)
        imshow('motifs_images/sc.png','InitialMagnification','fit');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',0,'enable','off'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','off'); set(handles.S3_0,'enable','off')
        set(handles.s2_init,'enable','off');set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','off')
        set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
        set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
        
        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','off','String',0);
        set(handles.s3in_in,'Enable','off','String',0);
        set(handles.ki2_in,'String',0);
        set(handles.ks1_in,'String',0.3);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','on','String',35);
        set(handles.y2_in,'Enable','on','String',0.06);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',2.5e-5);
        set(handles.km3_in,'Enable','off','String',0);
        set(handles.y3_in,'Enable','off','String',0);
        set(handles.kdec3_in,'Enable','off','String',0);
        set(handles.ks3_in,'Enable','off','String',0);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','off','String',0);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
        
        set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
        set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
        set(handles.gam_text1,'Enable','off');set(handles.gam_text2,'Enable','off')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1 - f_2X_2$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
        eqtx3='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        
        handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2');
        handles.var_names=strvcat('S1','X1','X2');
        xpos=0.7;
%         switch growth
%             case {'Hoh'}
%                 xpos=0.82;
%                 strB=[{'S1'},{'X1'},{'X2'}];
%                 for k=1:length(handles.steqs)
%                     strB=[strB,handles.steqs{k}];
%                 end
%                 
%                 [strC,indC]=unique(strB,'first');
%                 handles.var_names=char(strB(sort(indC)));
% 
%         end
        
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx5T);
        cla(handles.txax)
        axes(handles.txax)
        
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        text(xpos,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'Predation'
        set(handles.model_label,'String','5 ODE: Food Chain with Waste Product Inhibition');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','on'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='fci';
        axes(handles.motif_image)
        imshow('motifs_images/fcpi.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',1,'enable','on'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','on'); set(handles.S3_0,'enable','on')
        set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','on')
        set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
        set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');

        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','off','String',0);
        set(handles.s3in_in,'Enable','off','String',0);
        set(handles.ki2_in,'String',0.0000035);
        set(handles.ks1_in,'String',0.3);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0.1);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','on','String',35);
        set(handles.y2_in,'Enable','on','String',0.06);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',2.5e-5);
        set(handles.km3_in,'Enable','off','String',0);
        set(handles.y3_in,'Enable','off','String',0);
        set(handles.kdec3_in,'Enable','off','String',0);
        set(handles.ks3_in,'Enable','off','String',0);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','on','String',0.43);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
        set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
        set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
        set(handles.gam_text1,'Enable','off');set(handles.gam_text2,'Enable','off')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1I_3$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1I_3 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma_0(1-Y_1)f_1X_1I_3 - f_2X_2$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
        eqtx5='$\frac{dS_3}{dt} = -DS_3 + f_2X_2$';
        I3tx = '$I_3 = \frac{1}{1+\frac{S_3}{K_{i,3}}}$';
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        
        handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2','KI3');
        handles.var_names=strvcat('S1','X1','S2','X2','S3');
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5);
        cla(handles.txax)
        axes(handles.txax)
        
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        text(0.7,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'No_interaction'
        set(handles.model_label,'String','4 ODE: No Common Metabolites');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','on');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='ni';
        axes(handles.motif_image)
        imshow('motifs_images/ni.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','off'); set(handles.S3_0,'enable','off')
        set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','on')
        set(handles.s2in_in,'enable','on'); set(handles.text55,'enable','on')
        set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');

        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.ks1_in,'String',0.3);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','on','String',0);
        set(handles.s3in_in,'Enable','off','String',0);
        set(handles.ki2_in,'String',0);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','on','String',20);
        set(handles.y2_in,'Enable','on','String',0.06);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',0.5);
        set(handles.km3_in,'Enable','off','String',0);
        set(handles.y3_in,'Enable','off','String',0);
        set(handles.kdec3_in,'Enable','off','String',0);
        set(handles.ks3_in,'Enable','off','String',0);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','on','String',0.43);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
        
        set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
        set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
        set(handles.gam_text1,'Enable','off');set(handles.gam_text2,'Enable','off')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = D(S_{2,in} - S_2) - f_2X_2$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        
        handles.params_sim=strvcat('S1in','D','S2in','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2');
        handles.var_names=strvcat('S1','X1','S2','X2');
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
        cla(handles.txax)
        axes(handles.txax)
        
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        text(0.7,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'Cooperation'
        set(handles.model_label,'String','4 ODE: Syntrophy');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','on'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='syn';
        axes(handles.motif_image)
        imshow('motifs_images/syn.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','off'); set(handles.S3_0,'enable','off')
        set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','on')
        set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
        set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');

        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','off','String',0);
        set(handles.s3in_in,'Enable','off','String',0);
        set(handles.ki2_in,'String',0.0000035);
        set(handles.ks1_in,'String',0.3);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','on','String',35);
        set(handles.y2_in,'Enable','on','String',0.06);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',2.5e-5);
        set(handles.km3_in,'Enable','off','String',0);
        set(handles.y3_in,'Enable','off','String',0);
        set(handles.kdec3_in,'Enable','off','String',0);
        set(handles.ks3_in,'Enable','off','String',0);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','on','String',0.43);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
        set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
        set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
        set(handles.gam_text1,'Enable','off');set(handles.gam_text2,'Enable','off')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1I_2$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1I_2 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma_0(1-Y_1)f_1X_1I_2 - f_2X_2$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
        I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
        
        growth=handles.growthmodel;
        
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        if ~isempty(handles.thermeqs) && handles.thermeqs{1}=='Null'
            return
        end
      
        handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2','KI2');

        handles.var_names=strvcat('S1','X1','S2','X2');
        switch growth
            case {'Thermodynamic'}
                handles.var_names=strvcat('S1','X1','S2','X2',char(handles.steqs));
        end
        
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5T);
        cla(handles.txax)
        axes(handles.txax)
        if size(eqtext,1)>6
            fs=10;
        elseif size(eqtext,1)==6
            fs=12;
        else
            fs=14;
        end
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',fs)
        text(0.7,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'Amensalism'
        set(handles.model_label,'String','5 ODE: Waste Product Inhibition');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','on'); set(handles.threesp,'Checked','off');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='pi';
        axes(handles.motif_image)
        imshow('motifs_images/wpi.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',1,'enable','on'); set(handles.x3_check,'value',0,'enable','off')
        set(handles.s3_init,'enable','on'); set(handles.S3_0,'enable','on')
        set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','off')
        set(handles.X3_0,'enable','off'); set(handles.S2_0,'enable','on')
        set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
        set(handles.s3in_in,'enable','on'); set(handles.text48,'enable','on')
        set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');

        handles.plotsolution='time'; handles.plotdim='two';
        
        %Parameters
        set(handles.km1_in,'String',13);
        set(handles.y1_in,'String',0.04);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','off','String',0);
        set(handles.s3in_in,'Enable','on','String',0);
        set(handles.ki2_in,'String',0.0000035);
        set(handles.ks1_in,'String',0.3);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0.1);
        set(handles.x3_init,'String',0);
        set(handles.km2_in,'Enable','off','String',0);
        set(handles.y2_in,'Enable','off','String',0);
        set(handles.kdec2_in,'Enable','off','String',0);
        set(handles.ks2_in,'Enable','off','String',0);
        set(handles.km3_in,'Enable','on','Value',21);
        set(handles.y3_in,'Enable','on','String',0.04);
        set(handles.kdec3_in,'Enable','on','String',0.02);
        set(handles.ks3_in,'Enable','on','String',1e-3);
        set(handles.ks32_in,'Enable','off','String',0);
        set(handles.gamma0,'Enable','on','String',0.43);
        set(handles.gamma1,'Enable','off','String',0);
        set(handles.gamma2,'Enable','off','String',0); 
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.text14,'Enable','off');set(handles.text15,'Enable','off')
        set(handles.text16,'Enable','off');set(handles.text19,'Enable','off')
        set(handles.text43,'Enable','on');set(handles.text44,'Enable','on')
        set(handles.text45,'Enable','on');set(handles.text47,'Enable','on')
        set(handles.gam_text1,'Enable','off'); set(handles.gam_text2,'Enable','off')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma_0(1-Y_1)f_1X_1$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_3f_3X_2I_2 - k_{dec,3}X_2$';
        eqtx5='$\frac{dS_3}{dt} = D(S_{3,in}-S_3) - f_3X_2I_2$';
        I3tx = '$I_4 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        
        handles.params_sim=strvcat('S1in','D','S3in','km1','Ks1','Y1','km3','Ks3','Y3','kdec1','kdec2','KI2');
        handles.var_names=strvcat('S1','X1','S2','X2','S3');
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5);
        cla(handles.txax)
        axes(handles.txax)
        
        text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        text(0.7,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
        handles.eqtx=eqtext;
        
    case 'Threespecies'
        set(handles.model_label,'String','6 ODE: 3 species food-web');
        set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
        set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
        set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','on');
        set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
        handles.motif='ths';
        axes(handles.motif_image)
        imshow('motifs_images/ths.png');
        set(handles.plot_button,'enable','off')
        set(handles.plot_button2,'enable','off')
        set(handles.s1_check,'value',1,'enable','on'); set(handles.x1_check,'value',1,'enable','on')
        set(handles.s2_check,'value',1,'enable','on'); set(handles.x2_check,'value',1,'enable','on')
        set(handles.s3_check,'value',1,'enable','on'); set(handles.x3_check,'value',1,'enable','on')
        set(handles.s3_init,'enable','on'); set(handles.S3_0,'enable','on')
        set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','on')
        set(handles.X3_0,'enable','on'); set(handles.S2_0,'enable','on')
        set(handles.text55,'enable','on'); set(handles.text48,'enable','on')
        set(handles.text43,'enable','on'); set(handles.text44,'enable','on')
        set(handles.text45,'enable','on'); set(handles.text47,'enable','on')
        set(handles.text81,'enable','on'); set(handles.gam_text1,'enable','on')
        set(handles.gam_text2,'enable','on');
        
        %Parameters
        set(handles.km1_in,'String',29.12);
        set(handles.y1_in,'String',0.019);
        set(handles.kdec1_in,'String',0.02);
        set(handles.s1in_in,'String',5);
        set(handles.s2in_in,'Enable','on','String',0);
        set(handles.s3in_in,'Enable','on','String',0);
        set(handles.ki2_in,'String',0.0000035);
        set(handles.ks1_in,'String',0.000052);
        set(handles.d_in,'String',0.1);
        set(handles.time_in,'String',1000);
        set(handles.s1_init,'String',0.1);
        set(handles.x1_init,'String',0.1);
        set(handles.s2_init,'String',0.1);
        set(handles.x2_init,'String',0.1);
        set(handles.s3_init,'String',0.1);
        set(handles.x3_init,'String',0.1);
        set(handles.km2_in,'Enable','on','String',26);
        set(handles.y2_in,'Enable','on','String',0.04);
        set(handles.kdec2_in,'Enable','on','String',0.02);
        set(handles.ks2_in,'Enable','on','String',0.302);
        set(handles.km3_in,'Enable','on','Value',35);
        set(handles.y3_in,'Enable','on','String',0.06);
        set(handles.kdec3_in,'Enable','on','String',0.02);
        set(handles.ks3_in,'Enable','on','String',0.000025);
        set(handles.ks32_in,'enable','on','String',1e-6);
        set(handles.gamma0,'Enable','on','String',1.7069);
        set(handles.gamma1,'Enable','on','String',0.1429);
        set(handles.gamma2,'Enable','on','String',0.07069); 
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off')

        handles.plotsolution='time'; handles.plotdim='two';
        
        set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
        set(handles.text16,'Enable','on');set(handles.text19,'Enable','on'); set(handles.text55,'Enable','on')
        
        axes(handles.trajectoryplot)
        colorbar off
        %Equations
        eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
        eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
        eqtx3='$\frac{dS_2}{dt} = D(S_{2,in} - S_2) + \gamma_0(1-Y_1)f_1X_1 - f_2X_2I_2$';
        eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2I_2 - k_{dec,2}X_2$';
        eqtx5='$\frac{dS_3}{dt} = D(S_{3,in}-S_3) + \gamma_1(1-Y_2)f_2X_2I_2 - f_3X_3 - \gamma_2f_1X_1$';
        eqtx6='$\frac{dX_3}{dt} = -DX_3 + Y_3f_3X_3 - k_{dec,3}X_3$';
        I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';
        
        growth=handles.growthmodel;
        motif=handles.motif;
        dse=0; %Only set display equations
        growth_functions;
        if ~isempty(handles.thermeqs) && handles.thermeqs{1}=='Null'
            return
        end
        
        handles.params_sim=strvcat('S1in','D','S3in','S2in','km1','Ks1','Y1','km2','Ks2','Y2','km3','Ks3','Y3','kdec1','kdec2','kdec3','Ks3c','KI2');
        handles.var_names=strvcat('S1','X1','S2','X2','S3','X3',char(handles.steqs));
        eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5,eqtx6,eqtx5T);
        cla(handles.txax)
        axes(handles.txax)
        if size(eqtext,1)>6
            fs=10;
        else
            fs=12;
        end
        text(0,0.45,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',fs)
        text(0.81,0.5,handles.fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',12)
        handles.eqtx=eqtext;
end

motif_out=handles;



