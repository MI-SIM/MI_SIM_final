function [var]=set_variables(growth,gammas)
%set_variables - set the variables and parameters symbolically to be used in the model
%to analysis
%
% Author: Dr. Matthew Wade, School of Civil Engineering & Geosciences
% Newcastle University, Newcastle-upon-Tyne UK NE1 7RU
% email address: matthew.wade@ncl.ac.uk; dr.matthewwade@ncl.ac.uk
% alternative contact: Dr. Nick Parker, nick.parker@ncl.ac.uk
% Website: https://github.com/MI-SIM/MI-SIM
% September 2015; Last revision: 08-Mar-2016

syms S1 X1 S2 X2 S3 X3 km1 Ks1 Y1 kdec1 km2 Ks2 Y2 kdec2 km3 Ks3 Ks3c Y3 kdec3 KI2 D S1in S2in S3in Ks1a Ks2a Ks3a time1 gamma0 gamma1 gamma2 gamma3 gamma4 gamma5 gamma6 n1 n2 n3 S4 S5 S6 S7 S8

switch growth
    case {'Monod','Contois','Tessier','Haldane','Andrews'}
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2];
       
    case 'Moser'
        
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2 n1 n2 n3];
    
    case 'Thermodynamic'
        newgam={};
        for k=1:length(str2num(char(gammas)))
           newgam=[newgam,{['gamma',num2str(2+k)]}]; 
        end
        var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2 newgam]; 
        
%     case 'Hoh'
%         newgam={};
%         for k=1:length(str2num(char(gammas)))
%             newgam=[newgam,{['gamma',num2str(2+k)]}];
%         end
%         var=[km1 Y1 kdec1 km2 Y2 kdec2 km3 Y3 kdec3 KI2 D S1in S2in S3in time1 Ks1 Ks2 Ks3 Ks3c gamma0 gamma1 gamma2 newgam];
end

