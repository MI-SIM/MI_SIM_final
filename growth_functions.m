%growth_functions - generate the system growth functions
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIME
% February 2015; Last revision: 1-Apr-2016

%Define growth functions according to selection growth model

adt=1; f3=[]; eqtx5T=[]; %Only for thermodynamic
set(handles.gam_text0,'Enable','on')
switch growth
    case 'Monod'         %Monod growth function
        set(handles.gamma3,'Enable','off','String',0);
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}
                
                I2tx = '';
                if dse==1
                    f2=(km2*S2/(Ks2+S2));
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                f2tx = '$f_2 = \frac{k_{m,2}S_1}{K_{S,2} + S_1}$';
                if dse==1
                    f2=(km2*S1/(Ks2+S1));
                end
                
            case 'pi'
                f2tx = '';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
                if dse==1
                    f2=[];
                    f3=(km3*S3/(Ks3+S3));
                end
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                if dse==1
                    f2=(km2*S2/(Ks2+S2));
                end
                
            case 'ths'
                if dse==1
                    adt=S3/(Ks3c+S3);
                    f2=(km2*S2/(Ks2+S2));
                    f3=(km3*S3/(Ks3+S3));
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}\frac{S_3}{K_{S,3c} + S_3}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';
        end
        
        if dse==1
            f1=(km1*S1/(Ks1+S1))*adt;
        else
            set(handles.n1_in,'Enable','off','String',3)
            set(handles.n2_in,'Enable','off','String',3)
            set(handles.n3_in,'Enable','off','String',3)
            set(handles.text20,'Enable','off','String','n1')
            set(handles.text21,'Enable','off','String','n2')
            set(handles.text46,'Enable','off','String','n3')
        end
        
    case 'Contois'       %Contois growth function
        set(handles.gamma3,'Enable','off','String',0);
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{C,1}*X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{C,2}*X_2 + S_2}$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}
                
                I2tx = '';
                if dse==1
                    f2=(km2*S2/(Ks2*X2+S2));
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                f2tx = '$f_2 = \frac{k_{m,2}S_1}{K_{C,2}*X_2 + S_1}$';
                if dse==1
                    f2=(km2*S1/(Ks2*X2+S1));
                end
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                if dse==1
                    f2=(km2*S2/(Ks2*X2+S2));
                end
                
            case 'pi'
                f2tx = '';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{C,3}*X_3 + S_3}$';
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
                if dse==1
                    f2=[];
                    f3=(km3*S3/(Ks3*X3+S3));
                end
                
            case 'ths'
                if dse==1
                    adt=S3/(Ks3c*X3+S3);
                    f2=(km2*S2/(Ks2*X2+S2));
                    f3=(km3*S3/(Ks3*X3+S3));
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{C,1}*X_1 + S_1}\frac{S_3}{K_{C,3c}*X_3 + S_3}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{C,3}*X_3 + S_3}$';
                I2tx = '$I_3 = \frac{1}{1+\frac{S_3}{K_{i,3}}}$';
        end
        
        if dse==1
            set(handles.text46,'Enable','off','String','n3')
            f1=(km1*S1/(Ks1*X1+S1))*adt;
        else
            set(handles.n1_in,'Enable','off','String',3)
            set(handles.n2_in,'Enable','off','String',3)
            set(handles.n3_in,'Enable','off','String',3)
            set(handles.text20,'Enable','off','String','n1')
            set(handles.text21,'Enable','off','String','n2')
        end
        
    case 'Tessier'      %Tessier growth function
        set(handles.gamma3,'Enable','off','String',0);
        f1tx = '$f_1 = k_{m,1}\left(1-\exp\left[-\frac{S_1}{K_{S,1}}\right]\right)$';
        f2tx = '$f_2 = k_{m,2}\left(1-\exp\left[-\frac{S_2}{K_{S,2}}\right]\right)$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}
                
                I2tx = '';
                if dse==1
                    f2=km2*(1-exp(-S2/Ks2));
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                f2tx = '$f_2 = k_{m,2}\left(1-\exp\left[-\frac{S_1}{K_{S,2}}\right]\right)$';
                if dse==1
                    f2=km2*(1-exp(-S1/Ks2));
                end
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                if dse==1
                    f2=km2*(1-exp(-S2/Ks2));
                end
                
            case 'pi'
                f2tx = '';
                f3tx = '$f_3 = k_{m,3}\left(1-\exp\left[-\frac{S_3}{K_{S,3}}\right]\right)$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                if dse==1
                    f2=[];
                    f3=km3*(1-exp(-S3/Ks3));
                end
            case 'ths'
                if dse==1
                    adt=1-exp(-S3/Ks3c);
                    f2=km2*(1-exp(-S2/Ks2));
                    f3=km3*(1-exp(-S3/Ks3));
                end
                f1tx = '$f_1 = k_{m,1}\left(1-\exp\left[-\frac{S_1}{K_{S,1}}\right]\right)\left(1-\exp\left[-\frac{S_3}{K_{S,3c}}\right]\right)$';
                f3tx = '$f_3 = k_{m,3}\left(1-\exp\left[-\frac{S_3}{K_{S,3}}\right]\right)$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';
        end
        
        if dse==1
            f1=km1*(1-exp(-S1/Ks1))*adt;
        else
            set(handles.n1_in,'Enable','off','String',3)
            set(handles.n2_in,'Enable','off','String',3)
            set(handles.n3_in,'Enable','off','String',3)
            set(handles.text20,'Enable','off','String','n1')
            set(handles.text21,'Enable','off','String','n2')
            set(handles.text46,'Enable','off','String','n3')
        end
        
    case 'Moser'      %Moser growth function
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.n1_in,'Enable','on','String',3)
        set(handles.n2_in,'Enable','on','String',3)
        set(handles.n3_in,'Enable','off','String',3)
        set(handles.text20,'Enable','on','String','n1')
        set(handles.text21,'Enable','on','String','n2')
        set(handles.text46,'Enable','off','String','n3')
        
        f1tx = '$f_1 = k_{m,1}(\frac{S_1^{n_1}}{K_{S,1}+S_1^{n_1}})$';
        f2tx = '$f_2 = k_{m,2}(\frac{S_2^{n_2}}{K_{S,2}+S_2^{n_2}})$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}'
                
                I2tx = '';
                if dse==1
                    f2=km2*(S2^n2/(Ks2+S2^n2));
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                f2tx = '$f_2 = k_{m,2}(\frac{S_1^{n_2}}{K_{S,2}+S_1^{n_2}})$';
                if dse==1
                    f2=km2*(S1^n2/(Ks2+S1^n2));
                end
                
            case 'pi'
                f2tx = '';
                f3tx = '$f_3 = k_{m,3}(\frac{S_3^{n_3}}{K_{S,3}+S_3^{n_3}})$';
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                
                if dse==1
                    
                    f2=[];
                    f3=km3*(S3^n3/(Ks3+S3^n3));
                end
                
            case 'syn'
                
                I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';
                if dse==1
                    f2=km2*(S2^n2/(Ks2+S2^n2));
                end
                
            case 'ths'
                set(handles.n3_in,'Enable','on','String',3)
                set(handles.text46,'Enable','on','String','n3')
                if dse==1
                    adt=(S3^n3/(Ks3+S3^n3));
                    f2=km2*(S2^n2/(Ks2+S2^n2));
                    f3=km3*(S3^n3/(Ks3+S3^n3));
                end
                f1tx = '$f_1 = k_{m,1}(\frac{S_1^{n_1}}{K_{S,1}+S_1^{n_1}})(\frac{S_3^{n_3}}{K_{S,3c}+S_3^{n_3}})$';
                f3tx = '$f_3 = k_{m,3}(\frac{S_3^{n_3}}{K_{S,3}+S_3^{n_3}})$';
                I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';
        end
        
        if dse==1
            f1=km1*(S1^n1/(Ks1+S1^n1))*adt;
        end
        
    case 'Haldane'      %Haldane growth function
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.n1_in,'Enable','off','String',3)
        set(handles.n2_in,'Enable','off','String',3)
        set(handles.n3_in,'Enable','off','String',3)
        set(handles.text20,'Enable','off','String','n1')
        set(handles.text21,'Enable','off','String','n2')
        set(handles.text46,'Enable','off','String','n3')
        
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}
                
                I2tx = '';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1);
                    f2=km2*S2/(Ks2+S2);
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1);
                    f2=km2*S1/(Ks2+S1);
                end
                
            case {'pi'}
                
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1 + K_{i,2}S_1^2}$';
                f2tx = '';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                
                I2tx = '$I_2=1$';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1+(S1^2*KI2));
                    f2='';
                    f3=km3*S3/(Ks3+S3);
                end
            case {'syn'}
                
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1 + K_{i,2}S_1^2}$';
                f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
                
                I2tx = '$I_2=1$';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1+(S1^2*KI2));
                    f2=km2*S2/(Ks2+S2);
                end
                
            case 'ths'
                set(handles.n3_in,'Enable','on','String',3)
                set(handles.text46,'Enable','on','String','n3')
                if dse==1
                    adt=S3/(Ks3c+S3);
                    f1=(km1*S1/(Ks1+S1))*adt;
                    f2=km2*S2/(Ks2+S2+(S2^2*KI2));
                    f3=km3*S3/(Ks3+S3);
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}\frac{S_3}{K_{S,3c} + S_3}$';
                f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2 + K_{i,2}S_2^2}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2=1$';
        end
        
        
    case 'Andrews'      %Andrews growth function
        set(handles.gamma3,'Enable','off','String',0);
        set(handles.n1_in,'Enable','off','String',3)
        set(handles.n2_in,'Enable','off','String',3)
        set(handles.n3_in,'Enable','off','String',3)
        set(handles.text20,'Enable','off','String','n1')
        set(handles.text21,'Enable','off','String','n2')
        set(handles.text46,'Enable','off','String','n3')
        
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        f3tx = '';
        
        switch motif
            
            case {'fc','fci','ni'}
                
                I2tx = '';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1);
                    f2=km2*S2/(Ks2+S2);
                end
                
            case {'sc'}
                set(handles.gam_text0,'Enable','off')
                I2tx = '';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1);
                    f2=km2*S1/(Ks2+S1);
                end
                
            case {'pi'}
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1 + \frac{S_1^2}{K_{i,2}}}$';
                f2tx = '';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2=1$';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1+(S1^2/KI2));
                    f2=[];
                    f3=km3*S3/(Ks3+S3);
                end
                
            case {'syn'}
                
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1 + \frac{S_1^2}{K_{i,2}}}$';
                f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
                
                I2tx = '$I_2=1$';
                
                if dse==1
                    f1=km1*S1/(Ks1+S1+(S1^2/KI2));
                    f2=km2*S2/(Ks2+S2);
                end
                
            case 'ths'
                set(handles.n3_in,'Enable','on','String',3)
                set(handles.text46,'Enable','on','String','n3')
                if dse==1
                    adt=S3/(Ks3c+S3);
                    f1=(km1*S1/(Ks1+S1))*adt;
                    f2=km2*S2/(Ks2+S2+(S2^2/KI2));
                    f3=km3*S3/(Ks3+S3);
                end
                f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}\frac{S_3}{K_{S,3c} + S_3}$';
                f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2 + \frac{S_2^2}{K_{i,2}}}$';
                f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                I2tx = '$I_2=1$';
        end
        
    case 'Thermodynamic'    %Thermodynamic function
        out=[]; f3tx = '';
        switch motif
            
            case {'fc','sc','fci','ni','pi'}
                return
            case {'syn'}
                if dse~=1
                    gamma=str2num(get(handles.gamma0,'String')); exg=0;
                    while isempty(out)
                        [handles.thermeqs,eqtx5T,handles.gammas,handles.steqs,handles.dG,handles.dG_acc,handles.Temperature]=thermo_calc(handles.motif_name,gamma,growth,1);
                        if handles.thermeqs{1}=='Null'
                            return
                        end
                        out=1;
                    end
                    f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
                    f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
                    I2tx = ['$I_2=1-\exp\left({\frac{',handles.dG,'}{RT}}\right)$'];
                    %Set initial condition edit boxes
                    lst=length(handles.steqs);
                    if lst==2
                        set(handles.s3_init,'Enable','on','String',0.1)
                        set(handles.S3_0,'Enable','on','String',[handles.steqs{1},'(0)'])
                        set(handles.x3_init,'Enable','on','String',0.1)
                        set(handles.X3_0,'Enable','on','String',[handles.steqs{2},'(0)'])
                    elseif lst>2
                        set(handles.s3_init,'Enable','on','String',0.1)
                        set(handles.S3_0,'Enable','on','String',[handles.steqs{1},'(0)'])
                        set(handles.x3_init,'Enable','on','String',0.1)
                        set(handles.X3_0,'Enable','on','String',[handles.steqs{2},'(0)'])
                        for k=4:4+lst-3
                            eval(['set(handles.s',num2str(k),'_init,''Enable'',''on'',''String'',0.1)'])
                            eval(['set(handles.S',num2str(k),'_0,''Enable'',''on'',''String'',[handles.steqs{',num2str(k-1),'},''(0)''])'])
                        end
                    end
                    
                end
                
            case {'ths'}
                
                if dse==1
                    adt=S3/(Ks3c+S3);
                    f3=(km3*S3/(Ks3+S3));
                else
                    gamma=str2num(get(handles.gamma1,'String')); exg=1;
                    while isempty(out)
                        [handles.thermeqs,eqtx5T,handles.gammas,handles.steqs,handles.dG,handles.dG_acc,handles.Temperature]=thermo_calc(handles.motif_name,gamma,growth,1);
                        if handles.thermeqs{1}=='Null'
                            return
                        end
                        out=1;
                    end
                    f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}\frac{S_3}{K_{S,3c} + S_3}$';
                    f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
                    f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
                    I2tx = ['$I_2=1-\exp\left({\frac{',handles.dG,'}{RT}}\right)$'];
                    %Set initial condition edit boxes
                    lst=length(handles.steqs)-1;
                    for k=4:4+lst
                        eval(['set(handles.s',num2str(k),'_init,''Enable'',''on'',''String'',0.1)'])
                        eval(['set(handles.S',num2str(k),'_0,''Enable'',''on'',''String'',[handles.steqs{',num2str(k-3),'},''(0)''])'])
                    end
                end
                
        end
        
        if dse==1
            f1=(km1*S1/(Ks1+S1))*adt;
            f2=(km2*S2/(Ks2+S2));
        else
            set(handles.n1_in,'Enable','off','String',3)
            set(handles.n2_in,'Enable','off','String',3)
            set(handles.n3_in,'Enable','off','String',3)
            set(handles.text20,'Enable','off','String','n1')
            set(handles.text21,'Enable','off','String','n2')
            set(handles.text46,'Enable','off','String','n3')
            
            for nk=1:length(gamma)
                eval(['set(handles.gamma',num2str(nk-1+exg),',''Visible'',''on'',''enable'',''on'',''String'',handles.gammas(',num2str(nk),'))']);
            end
            
            for k=1:length(handles.gammas)-length(gamma)
                eval(['set(handles.gamma',num2str(k+2),',''Visible'',''on'',''enable'',''on'',''String'',handles.gammas(',num2str(k+length(gamma)),'))']);
                eval(['set(handles.gam_text',num2str(k+2),',''Visible'',''on'',''enable'',''on'')']);
            end
            set(handles.s3_init,'Enable','on','String',0.1)
            set(handles.S3_0,'Enable','on')
        end
        
%     case 'Hoh' %Hoh Model
%         f3tx = '';
%         
%         if dse==1
%             warning off
%             IG1a = (1-exp(sym(handles.dG{1}))); IG1b = (1+handles.kr{1}*exp(sym(handles.dG{1})));
%             IG2a = (1-exp(sym(handles.dG{2}))); IG2b = (1+handles.kr{2}*exp(sym(handles.dG{2})));
%             f1=km1*(S1*IG1a)/(Ks1+(S1*IG1b));
%             f2=km2*(S1*IG2a)/(Ks2+(S1*IG2b));
%         else
%             gamma=str2num(get(handles.gamma0,'String')); 
%             for xnum=1:2
%                 out=[];
%                 while isempty(out)
%                     if xnum==1
%                         [handles.thermeqs{xnum},eqtx,handles.gammas{xnum},handles.steqs{xnum},handles.dG{xnum},handles.dG_acc,handles.Temperature,handles.kr{xnum},handles.cmp_names,handles.gamma_count]=thermo_calc(handles.motif_name,gamma,growth,xnum);
%                     else
%                         [handles.thermeqs,eqtx,handles.gammas{xnum},handles.steqs{xnum},handles.dG{xnum},handles.dG_acc,handles.Temperature,handles.kr{xnum}]=thermo_calc(handles.motif_name,gamma,growth,xnum,handles.cmp_names,handles.thermeqs,eqtx,handles.gamma_count);
%                     end
%                     try
%                         if handles.thermeqs{1}=='Null'
%                             return
%                         end
%                     catch
%                         try
%                             if handles.thermeqs(1)=='Null'
%                                 return
%                             end
%                         catch
%                             if handles.thermeqs{1}{1}=='Null'
%                                 return
%                             end
%                             
%                         end
%                     end
%                     out=1;
%                 end
%                 
%             end
%             eqtx5T=eqtx;
%             I2tx = '';
%             
%             set(handles.n1_in,'Enable','off','String',3)
%             set(handles.n2_in,'Enable','off','String',3)
%             set(handles.n3_in,'Enable','off','String',3)
%             set(handles.text20,'Enable','off','String','n1')
%             set(handles.text21,'Enable','off','String','n2')
%             set(handles.text46,'Enable','off','String','n3')
%             
%             f1tx = '$f_1 = \frac{k_{m,1}S_1\Gamma_1^\alpha}{K_{S,1} + S_1\Gamma_1^\beta}$';
%             f2tx = '$f_2 = \frac{k_{m,2}S_1\Gamma_2^\beta}{K_{S,2} + S_1\Gamma_2^\beta}$';
%             I2tx = strvcat('$\Gamma_n^\alpha = (1-exp(\Delta G_n))$','$\Gamma_n^\beta = (1+k_{r,n}exp(\Delta G_n))$');
%             
%             set(handles.gam_text0,'Enable','off')
%             %Get number of gammas
%             len_gam=sum(cellfun('length',handles.gammas));
%             gams=[];
%             for kl = 1: length(handles.gammas)
%                 gams=[gams,cellfun(@str2num,handles.gammas{kl})];
%             end
%             for k=1:len_gam
%                 eval(['set(handles.gamma',num2str(k),',''Visible'',''on'',''enable'',''on'',''String'',gams(',num2str(k),'))']);
%                 eval(['set(handles.gam_text',num2str(k),',''Visible'',''on'',''enable'',''on'')']);
%             end
%             
%             %Set initial conditions dependent on variables
%             unqVars=unique([handles.steqs{1},handles.steqs{2}]);
%             for k=1:length(unqVars)
%                 eval(['set(handles.',lower(unqVars{k}),'_init,''Enable'',''on'',''String'',0.1)'])
%                 eval(['set(handles.',unqVars{k},'_0,''Enable'',''on'')'])
%             end
%             
%         end
        
end

if dse~=1
    handles.fnctext = strvcat(f1tx,f2tx,f3tx,I2tx);
end
