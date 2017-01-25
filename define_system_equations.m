function [eqs, f1, f2, f3, I2, I3, I4]=define_system_equations(motif,growth,handles)
%define_system_equations - generate the system of ODEs symbolically prior
%to analysis
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 22-Mar-2016

%Define symbolic variables (S1...X3) and parameters
syms S1 X1 S2 X2 km1 Ks1 Y1 kdec1 km2 Ks2 Y2 kdec2 KI2 D S1in S2in Ks1a Ks2a S3 X3 Y3 kdec3 Ks3a km3 Ks3 Ks3c S3in gamma0 gamma1 gamma2 gamma3 gamma4 gamma5 gamma6 n1 n2 n3 S4 S5 S6 S7 S8 dG

dse=1;

growth_functions;
R=8.3144598/1000;
try
    T=handles.Temperature;
end
%Define the equations of the system
switch motif
    
    case 'fc'
        
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        I2=[]; I3=[]; I4=[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        
    case 'sc'
        
        eq1=D*(S1in-S1)-f1*X1-f2*X2;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*X2+Y2*f2*X2-kdec2*X2;
        eq4=[]; eq5=[]; eq6=[];
        I2=[]; I3=[]; I4=[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        
        switch growth
            
            case 'Hoh'
                eqs=[eq1,eq2,eq3];
                for k=1:length(handles.thermeqs)
                    eval(['eq',num2str(3+k),'=eval(handles.thermeqs(',num2str(k),'));'])
                    eqs=[eqs,eval(['eq',num2str(3+k)])];
                end
                
        end
           
    case 'fci'
        
        I3=1/(1+S3/KI2);
        eq1=D*(S1in-S1)-f1*X1*I3;
        eq2=-D*X1+Y1*f1*X1*I3-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I3-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=-D*S3+f2*X2;
        eq6=[];
        I2=[]; I4=[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        
    case 'ni'
        
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        I2=[]; I3=[]; I4=[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        
    case 'syn'
        I2=1/(1+S2/KI2);
        switch growth
            case {'Haldane','Andrews'}
                I2=1;
                
            case 'Thermodynamic'
                I2=simplify(1-exp(sym(eval(handles.dG_acc/sym(R*T)))));
        end
        
        eq1=D*(S1in-S1)-f1*X1*I2;
        eq2=-D*X1+Y1*f1*X1*I2-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1*I2-f2*X2;
        eq4=-D*X2+Y2*f2*X2-kdec2*X2;
        eq5=[]; eq6=[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        switch growth

            case 'Thermodynamic'
                eqs=[eq1,eq2,eq3,eq4];
                for k=1:length(handles.thermeqs)
                    eval(['eq',num2str(4+k),'=eval(handles.thermeqs{',num2str(k),'});'])
                    eqs=[eqs,eval(['eq',num2str(4+k)])];
                end
                
        end
        I3=[]; I4=[];
        
    case 'pi'
        I4 = 1/(1+S2/KI2);
        switch growth
            case {'Haldane','Andrews'}
                I4=1;
        end
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=-D*S2+gamma0*(1-Y1)*f1*X1;
        eq4=-D*X2+Y3*f3*X2*I4-kdec2*X2;
        eq5=D*(S3in-S3)-f3*X2*I4;
        eq6=[];
        I2=[]; I3 =[];
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        
    case 'ths'
        I2=1/(1+S3/KI2);
        switch growth
            case {'Haldane','Andrews'}
                I2=1;
                
            case {'Thermodynamic'}
                I2=simplify(1-exp(sym(eval(handles.dG_acc/sym(R*T)))));
        end
        eq1=D*(S1in-S1)-f1*X1;
        eq2=-D*X1+Y1*f1*X1-kdec1*X1;
        eq3=D*(S2in-S2)+gamma0*(1-Y1)*f1*X1-f2*X2*I2;
        eq4=-D*X2+Y2*f2*X2*I2-kdec2*X2;
        eq5=D*(S3in-S3)+gamma1*(1-Y2)*f2*X2*I2-f3*X3-gamma2*f1*X1;
        eq6=-D*X3+Y3*f3*X3-kdec3*X3;
        eqs=[eq1,eq2,eq3,eq4,eq5,eq6];
        switch growth
            
            case 'Thermodynamic'
                for k=1:length(handles.thermeqs)
                    eval(['eq',num2str(5+k),'=eval(handles.thermeqs{',num2str(k),'});'])
                    eqs=[eqs,eval(['eq',num2str(5+k)])];
                end
        end
        
        I3=[]; I4=[];
        
end


