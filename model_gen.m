function [dP]=model_gen(time, init, parameters,h,h1,gfnc,growth,motif,flag,cno,thermeqs,gammas,dG,init_out,handles)
%model_gen - calls motif function to be solved with the ODE solver
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 24-Jan-2017

%Set precision
digits(32)

%Deal parameters
km1=parameters(1);
Y1=parameters(2);
kdec1=parameters(3);
km2=parameters(4);
Y2=parameters(5);
kdec2=parameters(6);
km3=parameters(7);
Y3=parameters(8);
kdec3=parameters(9);
KI2=parameters(10);
D=parameters(11);
S1in=parameters(12);
S2in=parameters(13);
S3in=parameters(14);
tt=parameters(15);
Ks1=parameters(16);
Ks2=parameters(17);
Ks3=parameters(18);
Ks3c=parameters(19);
gamma0=parameters(20);
gamma1=parameters(21);
gamma2=parameters(22);
switch growth
    case {'Moser'}
        n1=parameters(23);
        n2=parameters(24);
        n3=parameters(25);
    case {'Thermodynamic'}
        for k=1:length(gammas)
           eval(['gamma',num2str(2+k),'=gammas(',num2str(k),');']) 
        end
end

%Take the substrate and biomass values at each time point
outsiz=length(init);
S1=init(1);
X1=init(2);
S2=S1; 
I2=1; %Default inhibition term
I20=0;
dG0=0;
R=8.3144598/1000;
T=handles.Temperature;

%Set initial conditions dependent on motif chosen
try
    eval(init_out);
catch
    for k=1:outsiz
        eval([sprintf(handles.var_names(k,:)) '=init(k);'])
    end
end

%Evaluate growth function (symbolic -> numeric)

gnum=double(eval(gfnc));

%Generate system of ODEs for chosen motif
switch motif
    
    case 'fc'
        f1=gnum(1); f2=gnum(2);
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
    case 'sc'
        f1=gnum(1); f2=gnum(2);
        eq1=D*(S1in-S1)-f1*X1-f2*X2;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1; 
        eq3=-D*X2+Y2*f2*X2-kdec2*X2;
        eq4=[]; eq5=[]; eq6=[];
                
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
%         switch growth
%             case 'Hoh'
%                 dP=[eq1;eq2;eq3];
%                 for k=1:length(thermeqs)
%                     eval(['eq',num2str(3+k),'=thermeqs(',num2str(k),');'])
%                     try
%                         dP=[dP;eval(['eq',num2str(3+k)])];
%                     catch
%                         dP=[dP;['eq',num2str(3+k)]];
%                     end
%                 end
%                 dP=eval(dP);
%         end

    case 'fci'
        f1=gnum(1); f2=gnum(2);
        I3=gnum(3);
        eq1=D*(S1in-S1)-f1*X1*I3;
        eq2=-D*X1+Y1*f1*X1*I3-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I3-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=-D*S3+f2*X2; 
        eq6=[];
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
    case 'ni'
        f1=gnum(1); f2=gnum(2);
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
    case 'syn'
        f1=gnum(1); f2=gnum(2); I2=gnum(3);
        switch growth
            case 'Thermodynamic'
                try
                    I2 = 1-exp(eval(dG)/(R*T));
                catch
                    I2=1-exp(eval(dG/(R*T)));
                end
                dG0=eval(dG);
                I20=I2;
        end

        eq1=D*(S1in-S1)-f1*X1*I2;
        eq2=-D*X1+Y1*f1*X1*I2-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I2-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
        switch growth
            case 'Thermodynamic'
                dP=[eq1;eq2;eq3;eq4];
                for k=1:length(thermeqs)
                    eval(['eq',num2str(4+k),'=thermeqs{',num2str(k),'};'])
                    try
                        dP=[dP;eval(['eq',num2str(4+k)])];
                    catch
                        dP=[dP;['eq',num2str(4+k)]];
                    end
                end
                dP=eval(dP);
        end
        
    case 'pi'
        f1=gnum(1); f3=gnum(2);
        I4 = gnum(3);
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1;
        eq4=-D*X2+Y3*f3*X2*I4-kdec2*X2;
        eq5=D*(S3in-S3)-f3*X2*I4; 
        eq6=[];
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
    case 'ths'
        f1=gnum(1); f2=gnum(2); f3=gnum(3);
        I2=gnum(4); 
        
        switch growth
            case 'Thermodynamic'
                try
                    I2=1-exp(eval(dG/(R*T)));
                catch
                    I2 = 1-exp(eval(dG)/(R*T));
                end
                dG0=eval(dG);
                I20=I2;
        end
        
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)+gamma0*(1-Y1)*f1*X1-f2*X2*I2;
        eq4=-D*X2+Y2*f2*X2*I2-kdec2*X2;
        eq5=D*(S3in-S3)+gamma1*(1-Y2)*f2*X2*I2-f3*X3-gamma2*f1*X1; 
        eq6=-D*X3+Y3*f3*X3-kdec3*X3;
        dP=[eq1;eq2;eq3;eq4;eq5;eq6];
        switch growth
            case {'Thermodynamic'}
                dP=[eq1;eq2;eq3;eq4;eq5;eq6];
                for k=1:length(thermeqs)
                    eval(['eq',num2str(6+k),'=thermeqs(',num2str(k),');'])
                    try
                        dP=[dP;eval(eval(['eq',num2str(6+k),'{1}']))];
                    catch
                        try
                            dP=[dP;eval(['eq',num2str(6+k)])];
                        catch
                            dP=[dP;['eq',num2str(6+k)]];
                        end
                        dP=eval(dP);
                    end
                end
                
        end
end

%Update solver progress according to routine run

if isempty(flag)
    handles.timestamp=h;
    handles.progress=h1;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    total_tim=(time+1)*100/(tt+1);
    
    set(handles.timestamp,'String',['Time: ',TSMP])
    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
    
elseif flag==3
    handles.timestamp=h;
    handles.progress=h1;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    tm=cno+time;
    total_tim=(tm+1)*100/(tt+1);
    
    set(handles.timestamp,'String',['Time: ',TSMP])
    set(handles.progress,'String',['Progress: ',num2str(total_tim),'%'])
    
else
    handles.timestamp=h;
    c = clock; % [year month day hour minute seconds]
    TSMP = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
    set(handles.timestamp,'String',['Time: ',TSMP])
end
drawnow




