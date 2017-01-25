%stability_analysis - perform linear stability analysis and/or
%Routh-Hurwitz on the fixed-point solutions
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 11-Mar-2016

if exist('handles.htxt')
    delete(handles.htxt)
end
stabanaly=get(handles.lsanaly,'Value');
stabanaly2=get(handles.routhcrit,'Value');
jacobanaly=get(handles.jacobian_but,'Value');
set(handles.eigplot,'enable','off')
s={'-'}; s2={'-'};
%Get Jacobian of system
sym_vars = symvar(eqs);
var_new=[]; var_new_char=[];
for k=1:length(sym_vars)
    if strcmp(char(sym_vars(k)),'S1') || strcmp(char(sym_vars(k)),'X1') || strcmp(char(sym_vars(k)),'S2') || strcmp(char(sym_vars(k)),'X2') || strcmp(char(sym_vars(k)),'S3') || strcmp(char(sym_vars(k)),'X3')
        var_new=var_new;
        var_new_char=var_new_char;
    else
        var_new=[var_new;sym_vars(k)];
        var_new_char=strvcat(var_new_char,char(sym_vars(k)));
    end
end
for kk=1:length(var_new_char)
    eqs=subs(eqs,char(var_new(kk)),subs(var_new_char(kk,:)));
end
var_order=strvcat('S1','X1','S2','X2','S3','X3','S4','S5','S6');
kj=1; ind=[];
for k=1:length(sym_vars)
    if strmatch(char(sym_vars(k)),var_order)>0
        ind(kj)=strmatch(char(sym_vars(k)),var_order); kj=kj+1;
    end
end
%Create string elements
ind1=sort(ind);
var_order1=var_order(ind1,:);

%Create symbolic elements
vv=[];
for kh=1:length(ind1);
    eval(['v',num2str(kh),'= sym(var_order1(',num2str(kh),',:));']);
    eval(['vv=[vv,v',num2str(kh),'];']);
end

jac_sys=jacobian(vpa(eqs),vv);

%Create JPattern
[rws,cols]=size(jac_sys);
nuE=numel(jac_sys);

zons=ones(nuE,1);

zzer=find(jac_sys==0);

zons(zzer)=0;

Jpat=reshape(zons,rws,cols);

set(handles.spmatrix,'Visible','on');
axes(handles.spmatrix)
spy(Jpat,'r.',30)
axis tight
nz = nnz(Jpat);
pct = 100 / nuE;
xlabel(sprintf('nonzeros = %d (%.3f%%)',nz,nz*pct));
set(gca,'Xticklabel',cellstr(var_order1),'Yticklabel',cellstr(var_order1),'Xtick',1:length(Jpat),'Ytick',1:length(Jpat))
try
    delete(handles.htxt)
end
axes(handles.txax)

htxt=text(1.9,0.96,'Jacobian Sparsity Matrix','interpreter','latex','horiz','left','vert','middle','fontsize',12);
handles.htxt=htxt;
if stabanaly==1 %Linear stability analysis
    for i=1:length(fixed_numerical1(:,1))
        %Options for ODE solver (Jacobian not for ODE15 or ODE45)
        if jacobanaly==1
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[],'JPattern',Jpat);
        else
            options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        end
        
        %Get perturbation
        prt_s = get(handles.pert_p,'String');
        
        if isempty(str2num(prt_s))
            msgbox('The value entered for perturbation was not valid so the default value (1e-5) was used','Error','error')
            set(handles.pert_p,'String',0.00001);
            prt = 0.00001;
        else
            prt = str2num(prt_s);
        end
        
        %initial is the fixed point coordinates
        initial=['double([']; pertb=['['];
        
        for k=1:noeq
            initial=[initial,'fixed_numerical1(i,',num2str(k),') ']; 
            pertb=[pertb,'prt '];
        end
        
        initial=[initial,']);'];
        pertb=[pertb,'];'];
        
        initial=eval(initial);
        pertb=eval(pertb);
        
        init=initial+pertb;
        
        solver=char(strtrim(solver));
        set(handles.func_prog,'String',['Running: Stability Analysis ',num2str(i)],'ForegroundColor','r')
   
        h=handles.timestamp;
        h1=handles.progress;
        tt=time1;
        flag=[];cno=[];

        eval(['[~,yout2]=',solver,'(@model_gen, 0:stepsize:time1, init, options, parameters,h,h1,gfnc,growth,motif,flag,cno,handles.thermeqs,handles.gammas,handles.dG,init_out,handles);']);

        %calculate the difference between the trajectories
        dif1(:,i)=yout2(end,:)-initial;
        
        error=str2double(get(handles.error_val,'String'));
        if error<0 || isnumeric(error)==0
            set(handles.error_val,'String',1e-4);
            error=str2double(get(handles.error_val,'String'));
            msgbox('The value entered for Error was not valid so the default value was reset','Error','error')
        else
        end
        
        if norm(dif1(:,i))<error
            s(:,i)={'Stable'}; s2(:,i)={'-'};
        else
            s(:,i)={'Unstable'}; s2(:,i)={'-'};
        end
        set(handles.func_prog,'String',['Completed: Stability Analysis ',num2str(i)],'ForegroundColor',[0 0.6 1])
    end
end

if stabanaly2==1 %Routh-Hurwitz stability criterion
    
    for i=1:length(fixed_numerical1(:,1))
        jac_sysa=jac_sys;
        for kk=1:length(fixed_numerical1(i,:))
            var_n = fixed_numerical1(i,kk);
            if var_n==0
                var_n=1e-31;
            end
            jac_sysa=subs(jac_sysa,cellstr(var_order1(kk,:)),subs(var_n));
            
        end
        cpoly=charpoly(double(jac_sysa));
        eigen=roots(cpoly);
        if all(real(eigen)>0) && all(isreal(eigen))
            s2(:,i)={'Unstable node'};
        elseif all(real(eigen)<0) && all(isreal(eigen))
            s2(:,i)={'Stable node'};
        elseif any(real(eigen)>0) && any(real(eigen)<0) && all(isreal(eigen))
            s2(:,i)={'Unstable saddle point'};
        elseif all(~isreal(eigen)) && all(real(eigen)>0)
            s2(:,i)={'Unstable spiral'};
        elseif all(~isreal(eigen)) && all(real(eigen)<0)
            s2(:,i)={'Stable spiral'};
        elseif all(~isreal(eigen)) && any(real(eigen)==0)
            s2(:,i)={'Circle'};
        else
            s2(:,i)={'Unstable'};
        end
        
        qeig(i,:)=double(eigen);

        set(handles.eigplot,'enable','on')
        clear jac_sysa
        
    end
            handles.eigen=qeig;
end
guidata(hObject,handles)
