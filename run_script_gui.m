%run_script_gui - m-file for running the selected analysis routine [via the
%RUN button]
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 8-Apr-2016

%% Preamble
%Enable reporting
set(handles.report_menu,'Enable','on');
set(handles.reportpanel,'Visible','on'); set(handles.email_add,'Enable','on'); set(handles.server_add,'Enable','on')
set(handles.email_txt,'Enable','on'); set(handles.serv_txt,'Enable','on');

%Get value indicating if [jacobian sparsity matrix; normalised eigenvalues;
%plot overlay] will be used in the solver
jacobanaly=get(handles.jacobian_but,'Value');
ovra=get(handles.normeig,'Value');
ovro=get(handles.overlay,'Value');
gammas=[]; thermeqs=[]; dG=[]; init_out=[];
%For thermodynamics
switch handles.growthmodel
    case {'Thermodynamic'}
        thermeqs=handles.thermeqs;
        gammas=handles.gammas;
        dG=handles.dG_acc;
        
        gms=str2num(char(gammas));
        
        %Convert remaining gammas
        for k=1:length(gms)
            eval(['gamma',num2str(k+2),'=gms(',num2str(k),');']);
        end
        
end

%Clear figure axes and text information
if ovra==0
    try
        cla(handles.solutionplot)
        delete(handles.legend1)
    catch
        try
            for k=1:length(handles.hsub)
                cla(handles.hsub(k))
                delete(handles.hlegend(k))
            end
        end
    end
end

if ovro==0
    cla(handles.trajectoryplot)
end
if exist('handles.htxt')
    delete(handles.htxt)
end
if exist('handles.stbtxt')
    delete(handles.stbtxt)
end
%Define system equations
[eqs, f1, f2, f3, I2, I3, I4]=define_system_equations(motif,growth,handles);
gfnc=[f1,f2,f3,I2,I3,I4]; %Growth functions

%digits(32) %Set precision
format long
%% Algorithms
%Case structure for analysis routines
switch handles.simtype
    case 'single_p' %Single-point analysis
        
        %Get the symbolic variables and parameters
        variables=set_variables(growth,gammas);
        parameters=eval(variables);
        %Check number of equations
        noeq=length(eqs);
        
        %Substitute numerical values into equations
        set(handles.func_prog,'String','Running: Initialisation','ForegroundColor','r','Visible','on')
        set(handles.progress,'String',['Progress: ',num2str(0),'%'])
        drawnow
        for k=1:noeq
            eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
        end
        
        %% Fixed points solutions
        %Make assumptions (X,S real, non-negative)
        syms S1 X1 S2 X2 S3 X3 S4 S5 S6 S7 real
        assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); assumeAlso(X3>=0); assumeAlso(S4>=0); assumeAlso(S5>=0); assumeAlso(S6>=0); assumeAlso(S7>=0);
        
        %Use mupad to get all solutions
        eqnl=['['];
        clvr=cellstr(handles.var_names);
        eqvr=sym(clvr); 
        extra_init=['['];
        
        init_out=[];                               %Initial conditions initialisation

        for k=1:noeq
            eval(['eqnl=[eqnl,'' eq',num2str(k),'_numerical==0''];']);
            extra_init=[extra_init,[char(clvr(k)),'_init,']];
        end
        for k=3:noeq
            init_out=[init_out,[char(clvr(k)),'=init(',num2str(k),');']];            
        end
        eqnl=[eqnl,']'];
        extra_init(end)=']';
        dgts=str2num(get(handles.prec_num,'String'));
        digits(dgts)
        try
            sol=feval(symengine,'solve',simplify(eval(eqnl)),eqvr);
            
            try
                temp=vpa(sol);
                members=children(children(temp));
                solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
                num_sol=children(solutions{2});
                %Check for valid imaginary parts
                s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
                if s_imag==0
                    doub_sol=real(double(vpa(num_sol)));
                else
                    doub_sol=double(vpa(num_sol));
                end
            catch
                doub_sol=zeros(numel(sol),6);
                for i=1:numel(sol)
                    % Loop over the number of solutions
                    a=sol(i);
                    for j=1:numel(a)
                        % Loop over how many component
                        b=children(a(j));
                        doub_sol(i,j)=double(b(2));
                    end
                end
            end
        catch
            sol=feval(symengine,'solve',eval(eqnl),eqvr);
            try
                temp=vpa(sol);
                members=children(children(temp));
                solutions = cellfun(@(c)c(:),members,'UniformOutput',false);
                num_sol=children(solutions{2});
                %Check for valid imaginary parts
                s_imag=sum(any(imag(double(vpa(num_sol)))>1e-12));
                if s_imag==0
                    doub_sol=real(double(vpa(num_sol)));
                else
                    doub_sol=double(vpa(num_sol));
                end
            catch
                doub_sol=zeros(numel(sol),6);
                for i=1:numel(sol)
                    % Loop over the number of solutions
                    a=sol(i);
                    for j=1:numel(a)
                        % Loop over how many component
                        b=children(a(j));
                        doub_sol(i,j)=double(b(2));
                    end
                end
            end
        end
        digits 32
        clhead=clvr;
        
        init2=eval(extra_init);
        
        %Set ticks
        set(handles.s1_check,'Enable','off','Value',0); set(handles.x1_check,'Enable','off','Value',0); set(handles.s2_check,'Enable','off','Value',0);
        set(handles.x2_check,'Enable','off','Value',0); set(handles.s3_check,'Enable','off','Value',0); set(handles.x3_check,'Enable','off','Value',0);
        set(handles.s4_check,'Enable','off','Value',0); set(handles.s5_check,'Enable','off','Value',0); set(handles.s6_check,'Enable','off','Value',0);
        for kn=1:length(clhead)
            eval(['set(handles.',char(lower(clvr{kn})),'_check,''Enable'',''on'',''Value'',1,''String'',','''',clhead{kn},'''',')'])
        end

        %Remove invalid fixed points (fp<0)
        [ii,jj]=find(doub_sol<0);
        iu = unique(ii);
        doub_sol(iu,:)=[];
        
        [xx,yy]=size(doub_sol);
        if yy==1
            fixed_numerical1=doub_sol';
        else
            fixed_numerical1=doub_sol;
        end
        %Set the fixed points in the GUI
        lfn=length(fixed_numerical1(:,1));
        rnm={};
        for k=1:lfn
            rnm=[rnm,{['FP',num2str(k)]}];
        end
        set(handles.infotable,'Visible','on','Data',fixed_numerical1,'RowName',rnm,'ColumnName',clhead)
                
        %% Stability
        %Shows the stability of the fixed points in the GUI
        stability_analysis
        set(handles.stabtable,'Visible','on','Data',[s',s2'],'RowName',rnm)
        
        %Options for ODE solver
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'Jpattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        
        solver=char(strtrim(solver));
        
        %% Dynamics
        %Run solver - Dynamics
        set(handles.func_prog,'String','Running: Dynamics','ForegroundColor','r')
        h=handles.timestamp;
        h1=handles.progress;
        tt=time1; flag=[]; cno=[];
        te = cputime;
        eval(['[tout,yout]=',solver,'(@model_gen, [0:stepsize:time1], init2, options, parameters,h,h1,gfnc,growth,motif,flag,cno,handles.thermeqs,handles.gammas,handles.dG,init_out,handles);']);
        tfe= cputime-te;
        
        set(handles.func_prog,'String','Completed: Dynamics','ForegroundColor',[0 0.6 1])
        handles.tout=tout;
        handles.yout=yout;
        handles.fn=fixed_numerical1;
        
        switch growth
            
            case 'Thermodynamic'
                %Calculate dG and I2 for thermodynamics
                %Create Matlab function of dG
                symdG=sym(eval(handles.dG_acc));
                mfdG=matlabFunction(symdG);
                noS=length(symvar(symdG));
                yo='handles.dG0_out=mfdG(';
                
                indx=strmatch('S',clhead);
                
                for k=1:length(indx)
                    kk=indx(k);
                    yo=[yo,'yout(:,',num2str(kk),'),'];
                end
                yo(end:end+1)=');';
                eval(yo);
                handles.I20=1-exp(handles.dG0_out);
                fpstring=strvcat('Plots...','   Solutions','   deltaG','   Inhibition function');
                set(handles.furth_plots,'Visible','on','String',fpstring)
        end
        
        %Plot results
        which_plot='sol';
        plot_results;
        set(handles.plot_button,'enable','on')
        set(handles.plot_button2,'enable','on')
        set(handles.twodphase,'Enable','on'); set(handles.threedphase,'Enable','on')
        
    case 'multiple_p' %Multiple-point analysis
        %Plug in the values for numerical fixed points, ommitting selected
        %parameters (i.e. the parameter pair to sweep)
        
        variables = symvar(eqs);
        str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
        val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
        param1=str1(val1,:); param2=str2(val2,:);
        
        %Convert to cell
        var_cell=sym2cell(variables);
        str1=[]; sym_v=[];
        str_omit=[{'S1'},{'X1'},{'S2'},{'X2'},{'S3'},{'X3'},{'S4'},{'S5'},{'S6'},{'S7'},{strtrim(param1)},{strtrim(param2)}];
        ind=1;
        for k=1:length(var_cell)
            
            if sum(strcmp(char(var_cell{k}),str_omit))==0
                
                str1=strvcat(str1,char(var_cell{k}));
                pind(ind)=k;
                ind=ind+1;
            else
                str1=strvcat(str1,[]);
            end
        end
        variables=variables(pind);
        
        %Var_names
        for kk=1:length(handles.var_names)
            sym_v=[sym_v,sym(cellstr(handles.var_names(kk,:)))];
        end
        
        %Evaluate parameters to numeric
        for k=1:length(str1)
            parameters(k)=eval(str1(k,:));
        end
        
        noeq=length(eqs);
        
        for k=1:noeq
            eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
        end
        
        set(handles.progress,'String',['Progress: ',num2str(0),'%'])
        
        % Set fixed points empty
        S1_fix_num=[]; X1_fix_num=[]; S2_fix_num=[]; X2_fix_num=[]; S3_fix_num=[]; X3_fix_num=[]; S4_fix_num=[]; S5_fix_num=[]; S6_fix_num=[]; S7_fix_num=[];
        syms S1 X1 S2 X2 S3 X3 S4 S5 S6 S7 real
        assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); assumeAlso(X3>=0);  assumeAlso(S4>=0);  assumeAlso(S5>=0);
        assumeAlso(S6>=0); assumeAlso(S7>=0);
        eval(['syms ',param1]); eval(['syms ',param2]);
        r1=str2double(get(handles.min_p1,'String')); r2=str2double(get(handles.max_p1,'String')); s1a=str2double(get(handles.step_p1,'String'));
        r3=str2double(get(handles.min_p2,'String')); r4=str2double(get(handles.max_p2,'String')); s2a=str2double(get(handles.step_p2,'String'));
        range1=linspace(r1,r2,s1a);
        range2=linspace(r3,r4,s2a);
        
        eqs_n=subs(eqs,variables,parameters);
        jac_sys=jacobian(vpa(eqs_n),sym_v);
        
        xln=repmat(1:s1a,1,s2a);
        xlm=reshape(xln,s1a,s2a);
        ylm=repmat(1:s1a,s2a,1);
        yln=reshape(ylm,1,s1a*s2a);
        
        jacs=matlabFunction(jac_sys,'Vars',cellstr(strvcat(param1,param2,handles.var_names)));
        set(handles.func_prog,'String','Running: Multiple op. points - Mupad solver','ForegroundColor','r')
        
        eqtn=['['];
        for k=1:noeq
            eqtn=[eqtn,['eq',num2str(k),'_numerical==0 ']];
        end
        eqtn(end)=']';
        
        sc=cellstr(handles.var_names);
        
        drawnow
        sol=feval(symengine,'solve',vpa(eval(eqtn)),sc);
        eval(['syms ',param1]);
        eval(['syms ',param2]);
        
        drange=range2;
        srange=range1;
        [X,Y]=meshgrid(srange,drange);
        result_o = cell(s1a,s2a);
        
        for k=1:s2a
            
            temp = vpa( subs(sol,{param1,param2},{X(k,:),Y(k,:)}) );
            members = children(temp); %Find children of result
            solutions = cellfun(@(c)c(2),members,'UniformOutput',false); %Find numerical solutions
            r={};
            try
                result_o(k,:) = cellfun(@symToDouble,solutions,'UniformOutput',false);
            catch
                for kn=1:s2a
                    r{kn}={[X(k,kn);zeros(noeq-1,1)]};
                end
                result_o(k,:) = r;
            end
            
        end
        
        em_X=strfind(cellstr(handles.var_names),'X');
        no_X=sum(cell2mat(em_X));
        
        ff=find(~cellfun(@isempty,em_X));
        
        rsl=cellfun(@cell2mat,result_o,'UniformOutput',false);
        [rslw,rsll]=size(rsl);
        total_steps=s1a*s2a;
        no_steps=zeros(s1a,s2a);
        
        set(handles.func_prog,'String','Running: Multiple op. points - Calculating stability','ForegroundColor','r')
        drawnow
        for kk=1:rslw
            for ii=1:rsll
                ind=ii+s2a*(kk-1);
                total_tim=(ind)*100/(total_steps);
                set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
                c = clock; % [year month day hour minute seconds]
                TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
                set(handles.timestamp,'String',['Time: ',TSMP])
                drawnow
                Zr=rsl{kk,ii};
                nS=zeros(noeq,numel(Zr)/noeq);
                nS_a=[];
                for zk=1:noeq
                    nS_a=[nS_a,['nS(',num2str(zk),',i),']];
                end
                nS_a(end)='';
                for i=1:numel(Zr)/noeq
                    nS(:,i)=Zr(:,i);
                    jctxt=['jacs(range2(ylm(ii,kk)),',nS_a,');'];
                    jctxt2=['jacs(range1(xlm(ii,kk)),range2(ylm(ii,kk)),',nS_a,');'];
                    try
                        jac_P=eval(jctxt);
                    catch
                        jac_P=eval(jctxt2);
                    end
                    cpoly=charpoly(jac_P);
                    
                    eigen=roots(cpoly);
                    
                    if all(real(eigen<=0))==1
                        stab_k=1;
                    else
                        stab_k=0;
                    end
                    %Get existence and stability at point
                    fpn=[];
                    for k=1:no_X
                        fpn=[fpn,nS(ff(k),i)];
                    end
                    jk=0;
                    if nS(:,i)>=0
                        
                        if no_X==2
                            
                            if nnz(fpn)==2
                                jk=4;
                            elseif nnz(fpn)==0
                                jk=1;
                            elseif nnz(fpn)==1
                                nzu=find(fpn>0);
                                if nzu==1
                                    jk=2;
                                elseif nzu==2
                                    jk=3;
                                end
                            end
                            
                        elseif no_X==3
                            if nnz(fpn)==3
                                jk=6;
                            elseif nnz(fpn)==0
                                jk=1;
                            elseif nnz(fpn)==1
                                nzu=find(fpn>0);
                                if nzu==1
                                    jk=3;
                                elseif nzu==2
                                    jk=7;
                                elseif nzu==3
                                    jk=2;
                                end
                            elseif nnz(fpn)==2
                                nzu2=find(fpn>0);
                                if sum(nzu2)==3
                                    jk=4;
                                elseif sum(nzu2)==4
                                    jk=5;
                                elseif sum(nzu2)==5
                                    jk=8;
                                end
                            end
                            
                            %%%elseif no_X=4...etc
                            
                        end
                        
                        
                    end
                    
                    JmS = jk.*stab_k; jz=1-stab_k;
                    JmU = jk.*jz;
                    
                    if stab_k==1 && jk>0
                        Js(i,ind)=jk; Ss(i,ind) = {'S'}; Jn(i,ind)=jk;
                    elseif stab_k==0 && jk>0
                        Js(i,ind)=jk; Ss(i,ind) = {'U'}; Jn(i,ind)=jk+8.2;
                    else
                        Js(i,ind)=0; Ss(i,ind) = {'0'}; Jn(i,ind)=0;
                    end
                end
            end
        end
        
        sempS=cellfun(@isempty,Ss);
        Ss(sempS)={'0'};
        
        %Unique Js
        [unJs,unc,und]=unique(Js','rows');
        
        UniqueJs=Js(:,unc);
        
        %Unique Ss
        [unSs,uns,unt]=unique(cell2mat(Ss'),'rows');
        
        UniqueSs=Ss(:,uns);
        [m,n]=size(UniqueSs);
        %Find unique positions in SS
        
        out = sprintf('%.0f', Js) - '0';
        outU = unique(out);
        outF = outU(unique(outU)>0); ssF=[];
        for kF=1:length(outF)
            ssF=strvcat(ssF,['SS',num2str(outF(kF))]);
        end
        jjF=[];
        for k=1:n
            cc=find(unt==k);
            dd=Ss(:,cc);
            ee=Js(:,cc);
            [una,unb,unc]=unique(cell2mat(dd'),'rows');
            qF=cell.empty(0,length(outF));
            
            q=ee(:,1);
            q1=unique(q(q>0));
            q2=sort(q1);
            q_stab=[];
            for kk=1:length(q2)
                N=find(outF==q2(kk));
                q3=find(q==q1(kk));
                q_stab=dd(q3,1);
                if length(q3)>1
                    jstb=strjoin(q_stab);
                    q_stab=cellstr(strrep(jstb,' ','/'));
                end
                
                qF{N}=cell2mat(q_stab);
            end
            
            jjF=[jjF;strcat([['$\mathcal{J}_{',num2str(k),'}$'],qF])];
            
        end
        
        %%%%%MATCH JS AND SS%%%%%%%%%%
        
        sJ=sum(Jn,1);
        rJ=reshape(sJ,length(range1),length(range2));
        
        axes(handles.trajectoryplot)
        contourf(range1,range2,rJ','LineStyle','none','LevelStep',1)
        axis tight
        set(handles.colormp,'Visible','on');
        %LABELS
        xlabel(param1)
        ylabel(param2)
        
        %Populate stability matrix
        cla(handles.stabmatrix)
        axes(handles.stabmatrix)
        
        [a,b]=size(jjF);
        
        sm=cell(a+1,length(outF)+1);
        ssc=cellstr(ssF);
        ssct=ssc';
        sm(1,2:length(outF)+1)=ssct;
        sm(2:a+1,1:b)=jjF;
        [smr,smc]=size(sm);
        axes(handles.stabmatrix)
        [nm,nt]=size(sm);
        yn=ones(1,nt);
        xn=0.5;
        for j=1:nt
            xx(j)=(xn-0.16)+0.16;
            xn=xn+0.16;
        end
        for k=1:nm
            yy=yn*(0.8-(k-1)/10);
            if smr<9
                handles.stbtxt=text(xx+0.1,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',10);
            else
                handles.stbtxt=text(xx+0.1,yy,sm(k,:),'interpreter','latex','horiz','left','vert','middle','fontsize',6);
            end
        end
        
        %Change search arrays to numeric
        SsN=char(Ss); [rS,cS]=size(Ss);
        UsN=char(UniqueSs); [rU,cU]=size(UniqueSs);
        
        SsN(SsN=='S')='1'; UsN(UsN=='S')='1';
        SsN(SsN=='U')='2'; UsN(UsN=='U')='2';
        
        SsNu=str2num(SsN); UsNu=str2num(UsN);
        
        %Reshape and transpose
        SsNk=reshape(SsNu,rS,cS); UsNk=reshape(UsNu,rU,cU);
        SsNr=SsNk'; UsNr=UsNk';
        
        %Print locations of Jn
        xl=1:s1a;
        yl=1:s2a;
        [Xl Yl]=meshgrid(xl,yl); Xl=Xl'; Yl=Yl';
        axes(handles.trajectoryplot)
        for kj=1:cU
            
            Jna=ismember(SsNr,UsNr(kj,:),'rows');
            Jnj=reshape(Jna,s1a,s2a);
            
            jnj0 = mean(Xl(Jnj==1)); jnj1 = mean(Yl(Jnj==1));
            text(jnj0*range1(end)/length(range1),jnj1*range2(end)/length(range2),jjF(kj,1),'color','b','interpreter','latex','backgroundcolor','w','fontsize',14)
        end
        
        set(handles.func_prog,'String','Completed: Multiple operating points','ForegroundColor',[0 0.6 1])
        guidata(hObject,handles);
        
        
        %Save Figure to temp_fig
        fig_name=['temp_fig/',datestr(datetime),'_bifurcation_plot.pdf'];
        newfig_a=figure;
        axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
        hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
        close(newfig_a)
    case 'boa'
        basin_of_attraction;
        
    case 'pport'
        variables=set_variables(growth,handles.gammas);
        parameters=eval(variables);
        out_x_end=[]; out_y_end=[]; out_z_end=[];
        xlns=[]; ylns=[]; zlns=[]; out_x=[]; out_y=[]; out_z=[];
        
        %Check number of equations
        noeq=length(eqs);
        %For 3-dimensions
        is3D=get(handles.use_3v,'Value');
        
        %Options for ODE solver
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'Jpattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        
        solver=char(strtrim(solver));
        
        %%Get initial conditions
        min_ic1=str2num(get(handles.min_p1,'String')); min_ic2=str2num(get(handles.min_p2,'String'));
        max_ic1=str2num(get(handles.max_p1,'String')); max_ic2=str2num(get(handles.max_p2,'String'));
        step_ic1=str2num(get(handles.step_p1,'String'));
        
        if is3D==1
            min_ic3=str2num(get(handles.min_p3,'String')); max_ic3=str2num(get(handles.max_p3,'String'));
        end
        
        noeq=length(eqs);
        
        str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
        val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
        param1=str1(val1,:); param2=str2(val2,:);
        
        if is3D==1
            str3=get(handles.simparam3,'String'); val3=get(handles.simparam3,'Value');
            param3=str2(val3,:); init3(val3)=1;
        end
        
        init2=zeros(1,noeq);
        init2(val1)=1; init2(val2)=1;
        
        fz=find(init2==0);

        s1i=str2num(get(handles.s1_init,'String')); s2i=str2num(get(handles.s2_init,'String')); s3i=str2num(get(handles.s3_init,'String'));
        x1i=str2num(get(handles.x1_init,'String')); x2i=str2num(get(handles.x2_init,'String')); x3i=str2num(get(handles.x3_init,'String'));
        s4i=str2num(get(handles.s4_init,'String')); s5i=str2num(get(handles.s5_init,'String')); s6i=str2num(get(handles.s6_init,'String')); 

        vrs={'s1i','x1i','s2i','x2i','s3i','x3i','s4i','s5i','s6i'};
        
        for kj=1:noeq
            if any(fz==kj)
                init2(kj)=eval(char(vrs(kj)));
            end
        end
        
        %% Dynamics
        %Run solver - Dynamics at different initial conditions
        zt = length(0:stepsize:time1);cn=1;
        %Random or Fixed points
        grof=get(handles.step_p2,'Value');
        rng(0,'twister');
        if is3D==0
            tt=time1*step_ic1*step_ic1; flag=3; val3=[];
            parameters(find(variables=='time1'))=tt; %Change total time to correct value for loops
           % letr=param1(1);
            
%             %%%%%CHANGE THIS TO BE GENERIC%%%%%
%             Sarray=[1,3,5];
%             Xarray=[2,4,6];
%             if strcmp(letr,'S')
%                 s1m=find(Sarray==val1);
%                 s2m=find(Sarray==val2);
%                 Sarray([s1m,s2m])=[];
%                 val3=Sarray;
%             elseif strcmp(letr,'X')
%                 x1m=find(Xarray==val1);
%                 x2m=find(Xarray==val2);
%                 Xarray([x1m,x2m])=[];
%                 val3=Xarray;
%             end
%             
%             p3=str1(val3,:);
%             strp3=['handles.',lower(p3),'_init'];
%             init_val3=str2num(get(eval(strp3),'String'));
%             p3=str1(val3,:);
%             strp3=['handles.',lower(p3),'_init'];
%             init_val3=str2num(get(eval(strp3),'String'));

            if grof==1
                lsk1=linspace(min_ic1,max_ic1,step_ic1);
                lsk2=linspace(min_ic2,max_ic2,step_ic1);
                %lsk3=linspace(init_val3,init_val3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));
                %znl=repmat(lsk3,length(lsk3),1);
                %zlns=reshape(znl,1,numel(znl));
                
            elseif grof==2
                lsk1a = (max_ic1-min_ic1).*rand(1,step_ic1) + min_ic1; lsk1=sort(lsk1a);
                lsk2a = (max_ic2-min_ic2).*rand(1,step_ic1) + min_ic2; lsk2=sort(lsk2a);
               % lsk3=linspace(init_val3,init_val3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));
                %znl=repmat(lsk3,length(lsk3),1);
                %zlns=reshape(znl,1,numel(znl));
            end
            
        elseif is3D==1
            tt=time1*step_ic1*step_ic1*step_ic1; flag=3;
            if grof==1
                lsk1=linspace(min_ic1,max_ic1,step_ic1);
                lsk2=linspace(min_ic2,max_ic2,step_ic1);
                lsk3=linspace(min_ic3,max_ic3,step_ic1);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl)); zl=repmat(lsk3,1,length(lsk3)); zl2a=repmat(lsk3,length(lsk3),1);
                zl2b=reshape(zl2a,1,numel(zl2a)); zlns=meshgrid(zl,zl2b);
                
            elseif grof==2
                lsk1a = (max_ic1-min_ic1).*rand(1,step_ic1) + min_ic1; lsk1=sort(lsk1a);
                lsk2a = (max_ic2-min_ic2).*rand(1,step_ic1) + min_ic2; lsk2=sort(lsk2a);
                lsk3a = (max_ic3-min_ic3).*rand(1,step_ic1) + min_ic3; lsk3=sort(lsk3a);
                xlns=repmat(lsk1,1,length(lsk1)); ynl=repmat(lsk2,length(lsk2),1);
                ylns=reshape(ynl,1,numel(ynl));zl=repmat(lsk3,1,length(lsk3)); zl2a=repmat(lsk3,length(lsk3),1);
                zl2b=reshape(zl2a,1,numel(zl2a)); zlns=meshgrid(zl,zl2b);
            end
        end

        set(handles.func_prog,'String','Running: Dynamics','ForegroundColor','r')
        if is3D==0
            for k=lsk1
                for kk=lsk2
                    
                    init2(val1)=k;
                    init2(val2)=kk;
                    
                    h=handles.timestamp;
                    h1=handles.progress;
                    cno=(cn-1)*time1;
                    %eval(['[tout,yout]=',solver,'(@model_gen, [0:stepsize:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno,thermeqs,gammas,dG);']);
                    eval(['[tout,yout]=',solver,'(@model_gen, [0:stepsize:time1], init2, options, parameters,h,h1,gfnc,growth,motif,flag,cno,thermeqs,gammas,dG,init_out,handles);']);
                    out_x(1:zt,cn)=yout(:,val1);
                    out_y(1:zt,cn)=yout(:,val2);
                   % out_z(1:zt,cn)=yout(:,val3);
                    
                    out_x_end(cn)=yout(end,val1);
                    out_y_end(cn)=yout(end,val2);
                    %out_z_end(cn)=yout(end,val3);
                    
                    cn=cn+1;
                end
            end
            
        elseif is3D==1
            for k=lsk1
                for kk=lsk2
                    for kkk=lsk3
                        
                        init2(val1)=k;
                        init2(val2)=kk;
                        init2(val3)=kkk;
                        
                        h=handles.timestamp;
                        h1=handles.progress;
                        cno=(cn-1)*time1;
                        eval(['[tout,yout]=',solver,'(@model_gen, [0:stepsize:time1], init2, options, S1in, D, Y1, kdec1, Y2, kdec2, Y3, kdec3, km1, Ks1, km2, Ks2, km3, Ks3, Ks3c, KI2,gamma0, gamma1, gamma2, S2in, S3in,h,h1,tt,motif,flag,cno,thermeqs,gammas,dG);']);
                        
                        out_x(1:zt,cn)=yout(:,val1);
                        out_y(1:zt,cn)=yout(:,val2);
                        out_z(1:zt,cn)=yout(:,val3);
                        
                        out_x_end(cn)=yout(end,val1);
                        out_y_end(cn)=yout(end,val2);
                        out_z_end(cn)=yout(end,val3);
                        cn=cn+1;
                    end
                end
            end
        end
        
        set(handles.func_prog,'String','Completed: Dynamics','ForegroundColor',[0 0.6 1])
        axes(handles.trajectoryplot)
        
        %Separate branches
        rnd_out=round(out_x(end,:),5);
        un_rnd=unique(round(out_x(end,:),5));
        n=length(un_rnd);
        
        %Overlay plots?
        ovopt=get(handles.overlay,'Value');
        
        if n<7
            if ovopt==1
                hold on
                colr=jet(6);
                colre='y';
                mrk='-.';
            else
                hold off
                colr=lines(6);
                colre='g';
                mrk='-';
            end
        else
            if ovopt==1
                hold on
                colr=jet(n);
                colre='y';
                mrk='-.';
            else
                hold off
                colr=lines(n);
                colre='g';
                mrk='-';
            end
        end
        
        for k=1:n
            if is3D==0
                h=plot(out_x(:,rnd_out==un_rnd(k)),out_y(:,rnd_out==un_rnd(k)),'color',colr(k,:),'Linestyle',mrk);hold on
            elseif is3D==1
                h=plot3(out_x(:,rnd_out==un_rnd(k)),out_y(:,rnd_out==un_rnd(k)),out_z(:,rnd_out==un_rnd(k)),'color',colr(k,:),'Linestyle',mrk);hold on
            end
            set(h,'linewidth',2)
        end
        %Plot initial conditions
        if is3D==0
            plot(xlns,ylns,'ks','markersize',7,'markerfacecolor','r')
            %Plot steady-states (if available), otherwise plot final values
            plot(out_x_end,out_y_end,'ko','markersize',10,'markerfacecolor',colre)
            axis tight
            xlabel(param1); ylabel(param2)
            hold off
            box off
            
            %Save Figure to temp_fig
            fig_name=['temp_fig/',datestr(datetime),'_phaseportrait2d_plot.pdf'];
            newfig_a=figure;
            axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
            hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
            close(newfig_a)
        elseif is3D==1
            plot3(xlns,ylns,zlns,'ks','markersize',7,'markerfacecolor','r')
            %Plot steady-states (if available), otherwise plot final values
            plot3(out_x_end,out_y_end,out_z_end,'ko','markersize',10,'markerfacecolor',colre)
            axis tight
            xlabel(param1); ylabel(param2); zlabel(param3);
            hold off
            grid on
            box on
            set(handles.Dimensions,'Visible','off')
            
            %Save Figure to temp_fig
            fig_name=['temp_fig/',datestr(datetime),'_phaseportrait3d_plot.pdf'];
            newfig_a=figure;
            axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
            hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
            close(newfig_a)
        end
        
        set(handles.overlay,'Enable','on')
        handles.yout_phase=out_x; handles.yout_phase2=out_y; handles.yout_phase3=out_z;
        handles.yout_phase_end=out_x_end; handles.yout_phase_end2=out_y_end; handles.yout_phase_end3=out_z_end;
        handles.colr_phase=colr; handles.val3=val3; handles.mrk_phase=mrk;
        handles.colre_phase=colre;
        handles.xlns=xlns; handles.ylns=ylns; handles.zlns=zlns;
        set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
        set(gca,'uicontextmenu',handles.Dimensions)
        
end
handles.plothandle=get(gca,'Children');
%Get axes labels
handles.xlabelhandle=get(gca,'Xlabel');
handles.ylabelhandle=get(gca,'Ylabel');
handles.zlabelhandle=get(gca,'Zlabel');
guidata(hObject,handles)




