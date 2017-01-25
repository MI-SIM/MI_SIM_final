%basin_of_attraction - Plot basin of attraction
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 15-Mar-2016

%Get the values entered in the GUI

km1 = str2double(get(handles.km1_in,'String'));
Y1 = str2double(get(handles.y1_in,'String'));
kdec1 = str2double(get(handles.kdec1_in,'String'));
km2 = str2double(get(handles.km2_in,'String'));
Y2 = str2double(get(handles.y2_in,'String'));
kdec2 = str2double(get(handles.kdec2_in,'String'));
km3 = str2double(get(handles.km3_in,'String'));
Y3 = str2double(get(handles.y3_in,'String'));
kdec3 = str2double(get(handles.kdec3_in,'String'));
S1in = str2double(get(handles.s1in_in,'String'));
S2in = str2double(get(handles.s2in_in,'String'));
S3in = str2double(get(handles.s3in_in,'String'));
KI2 = str2double(get(handles.ki2_in,'String'));
Ks1 = str2double(get(handles.ks1_in,'String'));
Ks2 = str2double(get(handles.ks2_in,'String'));
Ks3 = str2double(get(handles.ks3_in,'String'));
Ks3c = str2double(get(handles.ks32_in,'String'));
n1 = str2double(get(handles.n1_in,'String'));
n2 = str2double(get(handles.n2_in,'String'));
n3 = str2double(get(handles.n3_in,'String'));
gamma0 = str2double(get(handles.gamma0,'String'));
gamma1 = str2double(get(handles.gamma1,'String'));
gamma2 = str2double(get(handles.gamma2,'String'));
gamma3 = str2double(get(handles.gamma3,'String'));
gamma4 = str2double(get(handles.gamma4,'String'));
gamma5 = str2double(get(handles.gamma5,'String'));
gamma6 = str2double(get(handles.gamma6,'String'));
D = str2double(get(handles.d_in,'String'));
time1 = str2double(get(handles.time_in,'String'));
S1_init=str2double(get(handles.s1_init,'String'));
X1_init=str2double(get(handles.x1_init,'String'));
S2_init=str2double(get(handles.s2_init,'String'));
X2_init=str2double(get(handles.x2_init,'String'));
S3_init=str2double(get(handles.s3_init,'String'));
X3_init=str2double(get(handles.x3_init,'String'));
S4_init=str2double(get(handles.s4_init,'String'));
S5_init=str2double(get(handles.s5_init,'String'));
S6_init=str2double(get(handles.s6_init,'String'));

%Get the values of the fixed points
[eqs, f1, f2, f3, I2, I3, I4]=define_system_equations(motif,growth,handles);

%Plug in the values for numeric fixed points
variables=set_variables(growth,handles.gammas);
parameters=eval(variables);

%Check number of equations
noeq=length(eqs);
eqst=[];
for k=1:noeq
    eval(['eq',num2str(k),'_numerical=subs(eqs(',num2str(k),'),variables,parameters);']);
    eqst=[eqst,['eq',num2str(k),'_numerical']];
end
%Solve system of ODEs
S1_fix_num=[]; X1_fix_num=[]; S2_fix_num=[]; X2_fix_num=[]; S3_fix_num=[]; X3_fix_num=[]; S4_fix_num=[]; S5_fix_num=[]; S6_fix_num=[]; S7_fix_num=[];
syms S1 X1 S2 X2 S3 X3 S4 S5 S6 S7 real
assumeAlso(S1>=0); assumeAlso(X1>=0); assumeAlso(S2>=0); assumeAlso(X2>=0); assumeAlso(S3>=0); assumeAlso(X3>=0);  assumeAlso(S4>=0);  assumeAlso(S5>=0);
assumeAlso(S6>=0); assumeAlso(S7>=0);

eqtn=['['];
for k=1:noeq
    eqtn=[eqtn,['eq',num2str(k),'_numerical==0 ']];
end
eqtn(end)=']';

sc=cellstr(handles.var_names);

sol=feval(symengine,'solve',eval(eqtn),sc);

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
    doub_sol=zeros(numel(sol),length(sc));
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

fixed_numerical1=doub_sol;

%Define parameter sweep
minp1a = get(handles.min_p1,'String');
minp1 = str2double(minp1a);
maxp1a = get(handles.max_p1,'String');
maxp1 = str2double(maxp1a);
stepp1a = get(handles.step_p1,'String');
stepp1 = str2double(stepp1a);
low_lim=minp1;
up_lim=maxp1;

%Vary the initial condition
x=linspace(low_lim,up_lim,stepp1);
y=linspace(low_lim,up_lim,stepp1);

%Label the axis
xvarls = get(handles.simparam1,'String'); yvarls=get(handles.simparam2,'String');
xvarop = get(handles.simparam1,'Value'); yvarop = get(handles.simparam2,'Value');

xvar = xvarls(xvarop,:); yvar = yvarls(yvarop,:);
axes(handles.trajectoryplot); cla

%Set up progress bar
total_steps = length(x)*length(y);
no_steps=zeros(length(x),length(y));

set(handles.func_prog,'String','Running: Basin of Attraction','ForegroundColor','r')

solver=char(strtrim(solver));
h=handles.timestamp; h1=[]; tt=[];

z_boa=zeros(length(x),length(y));
initial=zeros(length(x),length(y),2);  %set up a matrix of zeros
for i=1:length(x)
    for j=1:length(y)
        
        initial(i,j,1)=x(i);    %fill in the matrix with values
        initial(i,j,2)=y(j);    %fill in the matrix with values
        
        % varies two of the initial conditions and lets the user set the others

            if strcmp(xvar,'S1') && strcmp(yvar,'S1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S1')
                S1i=initial(i,j,1);
            elseif strcmp(yvar,'S1')
                S1i=initial(i,j,2);
            elseif sum(strcmp(sc,'S1'))==0
                S1i=[];
            else
                S1i=S1_init;
            end
            
            if strcmp(xvar,'X1') && strcmp(yvar,'X1')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X1')
                X1i=initial(i,j,1);
            elseif strcmp(yvar,'X1')
                X1i=initial(i,j,2);
            elseif sum(strcmp(sc,'X1'))==0
                X1i=[];
            else
                X1i=X1_init;
            end
            
            if strcmp(xvar,'S2') && strcmp(yvar,'S2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S2')
                S2i=initial(i,j,1);
            elseif strcmp(yvar,'S2')
                S2i=initial(i,j,2);
            elseif sum(strcmp(sc,'S2'))==0
                S2i=[];
            else
                S2i=S2_init;
            end
            
            if strcmp(xvar,'X2') && strcmp(yvar,'X2')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X2')
                X2i=initial(i,j,1);
            elseif strcmp(yvar,'X2')
                X2i=initial(i,j,2);
            elseif sum(strcmp(sc,'X2'))==0
                X2i=[];
            else
                X2i=X2_init;
            end
            
            if strcmp(xvar,'S3') && strcmp(yvar,'S3')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S3')
                S3i=initial(i,j,1);
            elseif strcmp(yvar,'S3')
                S3i=initial(i,j,2);
            elseif sum(strcmp(sc,'S3'))==0
                S3i=[];
            else
                S3i=S3_init;
            end
            
            if strcmp(xvar,'S4') && strcmp(yvar,'S4')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S4')
                S4i=initial(i,j,1);
            elseif strcmp(yvar,'S4')
                S4i=initial(i,j,2);
            elseif sum(strcmp(sc,'S4'))==0
                S4i=[];
            else
                S4i=S4_init;
            end
            
            if strcmp(xvar,'S5') && strcmp(yvar,'S5')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S5')
                S5i=initial(i,j,1);
            elseif strcmp(yvar,'S5')
                S5i=initial(i,j,2);
            elseif sum(strcmp(sc,'S5'))==0
                S5i=[];
            else
                S5i=S5_init;
            end
            
            if strcmp(xvar,'S6') && strcmp(yvar,'S6')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S6')
                S6i=initial(i,j,1);
            elseif strcmp(yvar,'S6')
                S6i=initial(i,j,2);
            elseif sum(strcmp(sc,'S6'))==0
                S6i=[];
            else
                S6i=S6_init;
            end
            
            if strcmp(xvar,'S7') && strcmp(yvar,'S7')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'S7')
                S7i=initial(i,j,1);
            elseif strcmp(yvar,'S7')
                S7i=initial(i,j,2);
            elseif sum(strcmp(sc,'S7'))==0
                S7i=[];
            else
                S7i=S7_init;
            end
            
            if strcmp(xvar,'X3') && strcmp(yvar,'X3')
                msgbox('Select two distinct variables','Error','error')
                return
                
            elseif strcmp(xvar,'X3')
                X3i=initial(i,j,1);
            elseif strcmp(yvar,'X3')
                X3i=initial(i,j,2);
            elseif sum(strcmp(sc,'X3'))==0
                X3i=[];
            else
                X3i=X3_init;
            end
     
            init1=[S1i, X1i, S2i, X2i, S3i, X3i, S4i, S5i, S6i, S7i];

        %Options for ODE solver
        options=odeset('RelTol',reltol,'AbsTol',abstol,'OutputFcn',[]);
        flag = 1;
        cno=[];
        handles.sc=sc;
        eval(['[tout,yout]=',solver,'(@model_gen, [0:stepsize:time1], init1, options, parameters,h,h1,gfnc,growth,motif,flag,cno,thermeqs,gammas,dG,init_out,handles);']);
        
        iter=j+length(y)*(i-1);
        total_tim=(iter)*100/(length(x)*length(y));
        set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
        drawnow
        
        %choose error
        error=str2double(get(handles.error_val,'String'));
        if error<0 || isnumeric(error)==0
            set(handles.error_val,'String',1e-6);
            error=str2double(get(handles.error_val,'String'));
            msgbox('The value entered for Error was not valid so the default value wa reset','Error','error')
        else
        end
        
        %Get basin of attraction
        for k=1:size(fixed_numerical1,1)
            if norm(yout(end,:)-fixed_numerical1(k,:))<error
                z_boa(i,j)=k-1;
            end
        end

    end
end
%plot the basin of attraction
axis([low_lim up_lim low_lim up_lim]);
colormap(hot); warning off
if mean(std(z_boa'))==0
else
    contourf(x,y,z_boa','linestyle','none','levelstep',1)
end

%Labels
xlabel(strcat(xvar,'(0)'),'FontSize',18)
ylabel(strcat(yvar,'(0)'),'FontSize',18)

hold on

for k=1:size(fixed_numerical1,1)
    %Label plot
    eval(['zons',num2str(k-1),'=zeros(length(x)*length(y),1);']);
    
    %Find indices
    eval(['zzer',num2str(k-1),'=find(z_boa==',num2str(k-1),');']);
    eval(['zons',num2str(k-1),'(zzer',num2str(k-1),')=1;']);
    
    %Reshape
    eval(['rM',num2str(k-1),'=reshape(zons',num2str(k-1),',length(x),length(y));']);
    
    xl=0:length(x)-1;
    yl=0:length(y)-1;
    
    [Xl Yl]=meshgrid(xl,yl);
    
    eval(['cX',num2str(k-1),'=mean(Xl(rM',num2str(k-1),'==1));']);
    eval(['cY',num2str(k-1),'=mean(Yl(rM',num2str(k-1),'==1));']);
    
    set(handles.func_prog,'String','Completed: Basin of Attraction','ForegroundColor',[0 0.6 1])
    
    eval(['ccX',num2str(k-1),'=cX',num2str(k-1),'*max(x)/length(x);']);
    eval(['ccY',num2str(k-1),'=cY',num2str(k-1),'*max(y)/length(y);']);
    
    eval(['text(ccX',num2str(k-1),',ccY',num2str(k-1),',''FP',num2str(k),''',''fontsize'',14,''color'',''w'',''backgroundcolor'',''k'',''fontweight'',''bold'');']);
    
end

%Save Figure to temp_fig
fig_name=['temp_fig/',datestr(datetime),'_boa_plot.pdf'];
newfig_a=figure;
axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
close(newfig_a)

