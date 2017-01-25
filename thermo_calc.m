function varargout = thermo_calc(varargin)
% THERMO_CALC MATLAB code for thermo_calc.fig
%      THERMO_CALC, by itself, creates a new THERMO_CALC or raises the existing
%      singleton*.
%
%      H = THERMO_CALC returns the handle to a new THERMO_CALC or the handle to
%      the existing singleton*.
%
%      THERMO_CALC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THERMO_CALC.M with the given input arguments.
%
%      THERMO_CALC('Property','Value',...) creates a new THERMO_CALC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before thermo_calc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to thermo_calc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help thermo_calc

% Last Modified by GUIDE v2.5 09-Oct-2016 14:12:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @thermo_calc_OpeningFcn, ...
                   'gui_OutputFcn',  @thermo_calc_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before thermo_calc is made visible.
function thermo_calc_OpeningFcn(hObject, eventdata, handles, varargin)

motif=varargin{1};
gammas=varargin{2};
set(handles.mass_text,'Enable','off') %Initial units - Molarity
handles.g_unit=2;
try
    growth=varargin{3};
    Xnum=varargin{4}; handles.xnum=Xnum;
    handles.growth=growth;
    handles.rte=Xnum-1;
    handles.gamma_count=1;
end

try
    handles.old_cmps=varargin{5};
    handles.old_eqs=varargin{6};
    handles.old_eqtx=varargin{7};
    handles.gamma_count=varargin{8};
end

if ~iscell(gammas)
    handles.gammas1=gammas;
    
    %Set gamma values
    for k=1:length(gammas)
        g=num2str(gammas(k));
        eval(['set(handles.p',num2str(k),'_gamma,''String'',g,''Enable'',''on'')']);
    end
    
else
    for k=1:length(gammas)
        eval(['set(handles.p',num2str(k),'_gamma,''Enable'',''on'')']);
    end
end
    
axes(handles.delta_G_ax)
text(0.13,0.65,'$\Delta G^0$ (kJ/mol)','interpreter','latex','horiz','left','vert','middle','fontsize',12)
axes(handles.delta_G_sel)
text(0.18,0.72,'$\Delta G_n^0$','interpreter','latex','horiz','left','vert','middle','fontsize',12)
axes(handles.axes5)
text(0.18,0.7,'$\Delta G^0$ (kJ/mol)','interpreter','latex','horiz','left','vert','middle','fontsize',12)
axes(handles.axes6)
text(0.5,0.74,'$\gamma_p$','interpreter','latex','horiz','left','vert','middle','fontsize',14)
switch motif
    
%     case 'Competition'
%         
%         str_react=['$$S1\mathop{\longrightarrow}^{X',num2str(Xnum),'}S2$$'];
%         set(handles.r1_text,'String','S1')
%         set(handles.p1_text,'String','S2')
%         handles.eqno=3; 
%         syms X1 Y1 f1 X2 Y2 f2
%         Xx=[X1,X2]; Ff=[f1,f2]; Yy=[Y1,Y2];
%         handles.funcN=(1-Yy(Xnum))*Ff(Xnum)*Xx(Xnum); 
%         handles.funcTa='+\gamma_';
%         handles.funcTb=['(1-Y_',num2str(Xnum),')f_',num2str(Xnum),'X_',num2str(Xnum)];
%         set(handles.react_sel,'Enable','off','Value',2)
%         set(handles.rr_text,'Visible','on'); set(handles.rr,'Visible','on');
%         set(handles.comp_title,'Visible','on')
%         handles.react_sy={'S1'};
%         handles.kr=1;
    
    case 'Cooperation'
        str_react='$$S1\mathop{\longrightarrow}^{X1}S2$$';
        set(handles.r1_text,'String','S1')
        set(handles.p1_text,'String','S2')
        handles.eqno=4;
        syms X1 Y1 I2 f1 
        handles.funcN=(1-Y1)*f1*X1*I2; 
        handles.funcTa='+\gamma_';
        handles.funcTb='(1-Y_1)f_1X_1I_2';
        %handles.kr=[];
        
    case 'Threespecies'
        str_react='$$S2\mathop{\longrightarrow}^{X2}S3$$';
        set(handles.r1_text,'String','S2')
        set(handles.p1_text,'String','S3')
        handles.eqno=6;
        
        syms X2 Y2 I2 f2 
        handles.funcN=(1-Y2)*f2*X2*I2; 
        handles.funcTa='+\gamma_';
        handles.funcTb='(1-Y_2)f_2X_2I_2';
        %handles.kr=[];
end
handles.motif=motif;
handles.cmp_names=[];
axes(handles.eq_disp)
text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
handles.output = hObject;
handles.systems={'H-O','H-O-N','H-O-S','H-O-N-S-C_inorganic','H-O-C_organic','H-O-S-C_metals/minerals','H-O-P','H-O-Halogens'};
handles.compounds1={'O2(g)','O2(aq)','H2O2(aq)','HO2(-)','H2O(l)','H2O(g)','OH(-)','H(+)','H2(g)','H2(aq)'};
handles.compounds2={'NO3(-)','HNO3(aq)','NO2(-)','HNO2(aq)','NO(g)','NO(aq)','N2O(g)','N2O(aq)','N2(g)','N2(aq)','NH3(g)','NH3(aq)','NH4(+)'};
handles.compounds3={'SO4(2-)','HSO4(-)','SO3(2-)','HSO3(-)','SO2(aq)','S2O3(2-)','HS2O3(-)','H2S2O3(aq)','S2O4(2-)','HS2O4(aq)','S2O5(2-)','S2O6(2-)','S2O8(2-)',...
    'S3O6(2-)','S4O6(2-)','S5O6(2-)','S(s)','HS(-)','H2S(aq)','H2S(g)','S2(g)','S2(2-)','S3(2-)','S4(2-)','S5(2-)'};
handles.compounds4={'CO2(g)','CO2(aq)','CO3(2-)','HCO3(-)','COS(g)','CO(g)','CO(aq)','CN(-)','HCN(aq)','OCN(-)','SCN(-)','CH4(g)','CH4(aq)'};
handles.compounds5={'Formate','Formic Acid(aq)','Acetate(-)','Acetic Acid(aq)','Glycolate(-)','Glycolic Acid(aq)','Propanoate(-)','Propanoic acid(aq)','Lactate(-)','Lactic Acid(aq)','Butanoic acid(aq)',...
    'Butanoate(-)','Pentanoic acid(aq)','Pentanoate(-)','Benzoate(-)','Benzoic Acid(aq)','Oxalate(2-)','H-Oxalate(-)','Oxalic Acid(aq)','Malonate(2-)','H-Malonate(-)','Malonic Acid(aq)',...
    'Succinate(2-)','H-Succinate(-)','Succinic Acid(aq)','Glutaric Acid(aq)','H-Glutarate(-)','Glutarate(2-)','Methanol(aq)','Ethanol(aq)','Propanol(aq)','2-Propanol(aq)','Butanol(aq)','Pentanol(aq)',...
    'Ethane(aq)','Propane(aq)','Butane(aq)','Pentane(aq)','Octane(l)','Nonane(l)','Decane(l)','Undecane(l)','Hexadecane(l)','Alanine(aq)','Arginine(aq)','Arginine(+)','Asparagine(aq)','Aspartic Acid(aq)',...
    'Aspartate(-)','Cysteine(aq)','Glutamic Acid(aq)','Glutamate(-)','Glutamine(aq)','Glycine(aq)','Histidine(aq)','Histidine(+)','Isoleucine(aq)','Leucine(aq)','Lysine(aq)','Lysine(+)','Methionine(aq)',...
    'Phenylalanine(aq)','Proline(aq)','Serine(aq)','Threonine(aq)','Tryptophan(aq)','Tyrosine(aq)','Valine(aq)','Methanamine(aq)','Toluene(aq)','Toluene(l)','Ethylbenzene(aq)','Ethylbenzene(l)'};
handles.compounds6={'Mg(2+)','MgOH(-)','MgHCO3(+)','MgCO3(aq)','Magnesite','Ca(2+)','CaOH(+)','CaHCO3(+)','CaCO3(aq)','Calcite','Dolomite','VO4(3-)','HVO4(2-)','H2VO4(-)','H3VO4(aq)','VO2(+)','VOOH(+)',...
    'VO(2+)','VO(+)','VOH(2+)','V(3+)','VOH(+)','V(2+)','Cr2O7(2-)','CrO4(2-)','HCrO4(-)','CrO(+)','CrOH(2+)','Cr(3+)','MnO4(-)','MnO4(2-)','Mn(3+)','Mn(2+)','MnOH(+)','MnO(aq)','MnO2(2-)','HMnO2(-)','Alabandite',...
    'Rhodochrosite','Fe(3+)','FeOH(2+)','FeO(+)','HFeO2(aq)','FeO2(-)','Hematite','Magnetite','Fe(2+)','FeOH(+)','FeO(aq)','HFeO2(-)','Pyyrhotite','Siderite','Pyrite','Co(3+)','CoOH(2+)','Co(2+)','CoOH(+)','CoO(aq)',...
    'HCoO2(-)','CoO2(2-)','Ni(2+)','NiOH(+)','NiO(aq)','HNiO2(-)','NiO2(2-)','Cu(2+)','CuOH(+)','CuO(aq)','HCuO2(-)','CuO2(2-)','Chalcopyrite','Covelite','Bornite','Copper','Zn(2+)','ZnOH(+)','ZnO(aq)','ZnO2(2-)',...
    'HZnO2(-)','Sphalerite','AsO4(3-)','HAsO4(2-)','H2AsO4(-)','H3AsO4(aq)','AsO2(-)','HAsO2(aq)','SeO4(2-)','HSeO4(-)','SeO3(2-)','HSeO3(-)','H2SeO3(aq)','Selenium','HSe(-)','MoO4(2-)','HMoO4(-)','Molybdenite','Ag(+)',...
    'Silver','WO4(2-)','HWO4(-)','Au(3+)','Au(+)','Gold','Hg(2+)','HgOH(+)','HgO(aq)','HHgO2(-)','Quicksilver','Pb(2+)','PbOH(+)','PbO(aq)','Galena','Anglesite','HPbO2(-)','UO2(2+)','UO2OH(+)','UO3(aq)','HUO4(-)','UO4(2-)',...
    'UO2(+)','UO2OH(aq)','UO3(-)','U(4+)','UOH(3+)','UO(2+)','HUO2(+)','UO2(aq)','Uraninite','HUO3(-)'};
handles.compounds7={'H3PO4(aq)','H2PO4(-)','HPO4(2-)','PO4(3-)','H4P2O7(aq)','H3P2O7(-)','H2P2O7(2-)','HP2O7(3-)','P2O7(4-)','H3PO2(aq)','H2PO2(-)','H3PO3(aq)','H2PO3(-)','HPO3(2-)'};
handles.compounds8={'HCl(aq)','Cl(-)','ClO(-)','HClO(aq)','ClO2(-)','HClO2(aq)','ClO3(-)','ClO4(-)','Br(-)','BrO(-)','HBrO(aq)','BrO3(-)','BrO4(-)','I(-)','IO(-)','HIO(aq)','IO3(-)','HIO3(aq)','IO4(-)'};
handles.temperatures=[2,18,25,37,45,55,70,85,100,115,150,200]+273.15; %Temperatures in Kelvin

set(handles.r1_cs,'Enable','off','String','Compounds')
set(handles.r1_ss,'Enable','off','String','Systems')
set(handles.r1_gas,'Value',0);
set(handles.r1_dg,'Visible','off','String','.');

set(handles.r2_cs,'Enable','off','String','Compounds')
set(handles.r2_ss,'Enable','off','String','Systems')
set(handles.r2_gas,'Value',0);
set(handles.r2_dg,'Visible','off','String','.');
set(handles.r2_text,'Enable','off');

set(handles.r3_cs,'Enable','off','String','Compounds')
set(handles.r3_ss,'Enable','off','String','Systems')
set(handles.r3_gas,'Value',0);
set(handles.r3_dg,'Visible','off','String','.');
set(handles.r3_text,'Enable','off');

set(handles.p1_cs,'Enable','off','String','Compounds')
set(handles.p1_ss,'Enable','off','String','Systems')
set(handles.p1_gas,'Value',0);
set(handles.p1_dg,'Visible','off','String','.');

set(handles.p2_cs,'Enable','off','String','Compounds')
set(handles.p2_ss,'Enable','off','String','Systems')
set(handles.p2_gas,'Value',0);
set(handles.p2_dg,'Visible','off','String','.');

set(handles.p3_cs,'Enable','off','String','Compounds')
set(handles.p3_ss,'Enable','off','String','Systems')
set(handles.p3_gas,'Value',0);
set(handles.p3_dg,'Visible','off','String','.');

set(handles.p4_cs,'Enable','off','String','Compounds')
set(handles.p4_ss,'Enable','off','String','Systems')
set(handles.p4_gas,'Value',0);
set(handles.p4_dg,'Visible','off','String','.');

set(handles.p5_cs,'Enable','off','String','Compounds')
set(handles.p5_ss,'Enable','off','String','Systems')
set(handles.p5_gas,'Value',0);
set(handles.p5_dg,'Visible','off','String','.');

set(handles.sel_comp_text,'Enable','off')
set(handles.sel_sys_text,'Enable','off')
set(handles.dG_ax,'Visible','off')
cla(handles.dG_ax)

set(handles.add_sn,'Enable','off')
set(handles.finish_but,'Enable','off')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes thermo_calc wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = thermo_calc_OutputFcn(hObject, eventdata, handles) 
try
    varargout{1}= handles.output1;
    varargout{2}= handles.output2;
    varargout{3}= handles.output3;
    varargout{4}= handles.output4;
    varargout{5}= handles.dG;
    varargout{6}= handles.dG_acc;
    varargout{7}= handles.Temperature;
    %varargout{8}= handles.kr;
    %varargout{9}= handles.cmp_names;
    %varargout{10}= handles.gamma_count;
catch
    varargout{1}={'Null'};
    varargout{2}=[];
    varargout{3}=[];
    varargout{4}=[];
    varargout{5}=[];
    varargout{6}=[];
    varargout{7}=[];
    %varargout{8}=[];
    %varargout{9}=[];
    %varargout{10}=[];
end
delete(gcf)

function temp_in_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function temp_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p1_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p1_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p2_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p2_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p3_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p4_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p4_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p5_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p5_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in product_select.
function product_select_Callback(hObject, eventdata, handles)
sel_num=get(handles.product_select,'Value');
rsel_num=get(handles.react_sel,'Value');
cla(handles.eq_disp)
if sel_num==1
    sel_num=2;
end

sel_numa=sel_num;
dd=2;

switch handles.motif

%     case {'Competition'}
%         rsym={'S1'};
%         psym={''};
%         num_start=3;
%         rsel_num=rsel_num-2; k=1; dd=1;
%         set(handles.cmp1,'Enable','on','Visible','on')
        
    case {'Cooperation'}
        rsym={'S1'};
        psym={'S2'};
        num_start=3;
        rsel_num=rsel_num-2; k=2;
        set(handles.cmp1,'Enable','on','Visible','on')
    case 'Threespecies'
        rsym={'S2'};
        psym={'S3'};
        rsel_num=rsel_num-2;
        if sel_num==1
            sel_num=2;
        end
        num_start=4; k=3;
end

if rsel_num<0
    rsel_numa=1
else
    rsel_numa=rsel_num;
end
for kk=2:sel_numa
    eval(['set(handles.p',num2str(kk),'_text,''Enable'',''on'',''String'',[''S'',num2str(num_start+kk-2+rsel_numa)])'])
    eval(['set(handles.p',num2str(kk),'_in,''Enable'',''on'')'])
    if handles.g_unit<3
        eval(['set(handles.p',num2str(kk),'_mass,''Enable'',''off'')'])
    else
        eval(['set(handles.p',num2str(kk),'_mass,''Enable'',''on'')'])
    end
    eval(['set(handles.p',num2str(kk),'_gas,''Enable'',''on'')'])
    eval(['set(handles.p',num2str(kk),'_gamma,''Enable'',''on'')'])
    eval(['set(handles.cmp',num2str(kk),',''Enable'',''on'',''Visible'',''on'')'])
end

if sel_numa<5
    for kk=sel_numa+1:5
        eval(['set(handles.p',num2str(kk),'_text,''Enable'',''off'')'])
        eval(['set(handles.p',num2str(kk),'_in,''Enable'',''off'')'])
        eval(['set(handles.p',num2str(kk),'_mass,''Enable'',''off'')'])
        eval(['set(handles.p',num2str(kk),'_gas,''Enable'',''off'')'])
        eval(['set(handles.p',num2str(kk),'_gamma,''Enable'',''off'')'])
        eval(['set(handles.cmp',num2str(kk),',''Enable'',''off'',''Visible'',''off'')'])
    end
end

rnsym={};
if rsel_num>0
    for kj=num_start:num_start+rsel_num-1
        rsym=[rsym,{['+S',num2str(kj)]}];
        rnsym=[rnsym,{['S',num2str(kj)]}];
    end
end

pnsym={};

z_s=k+1;
for kk=z_s:z_s+sel_num-dd
    psym=[psym,{['+S',num2str(kk)]}];
    pnsym=[pnsym,{['S',num2str(kk)]}];
end

handles.react_sy=rnsym;
handles.prod_sy=pnsym;
    switch handles.motif
        
%         case {'Competition'}
%             
%             str_react=['$$',cell2mat(rsym),'\mathop{\longrightarrow}^{X',num2str(handles.xnum),'} ',cell2mat(psym),'$$'];
%             axes(handles.eq_disp)
%             text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
          
        case {'Cooperation'}
            
            str_react=['$$',cell2mat(rsym),'\mathop{\longrightarrow}^{X1} ',cell2mat(psym),'$$'];
            axes(handles.eq_disp)
            text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
 
        case 'Threespecies'
            
            str_react=['$$',cell2mat(rsym),'\mathop{\longrightarrow}^{X2} ',cell2mat(psym),'$$'];
            axes(handles.eq_disp)
            text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
            
    end

handles.sel_num=sel_num;
guidata(hObject,handles)
    

% --- Executes during object creation, after setting all properties.
function product_select_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function r1_mass_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function r1_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p1_mass_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p1_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p2_mass_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p2_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p3_mass_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function p3_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p4_mass_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function p4_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function p5_mass_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function p5_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in calculate_but.
function calculate_but_Callback(hObject, eventdata, handles)
%Calculate Thermodynamic Inhibition

dg0=str2num(get(handles.dG0_in,'String'));
if isempty(dg0) | ~isnumeric(dg0)
    dg0=0;
end
R = 8.3144598/1000;                              %Universal gas constant
T = str2num(get(handles.temp_in,'String'));      %Temperature
handles.Temperature=T;
%Extract numerical value of substrate
b=get(handles.p1_text,'String');
c=sscanf(b,'%s%f');
Sd=str2num(char(c(2)));
Sdn=Sd+1;
Nn=1;

%For Reactants
mr1=str2num(get(handles.r1_in,'String'));
if handles.g_unit==3 %For COD conversion
    Mr1=str2num(get(handles.r1_mass,'String'));
    eqR=Mr1^mr1;
end

Symbr={get(handles.r1_text,'String')}; mr_store=mr1;
try
    dgR1=mr1*-str2num(get(handles.r1_dg,'String')); %Negative on product side
catch
    dgR1=1;
end
Rrs=get(handles.react_sel,'Value');
jg=2; rg_store=dgR1; 
if Rrs>2
    for kg=2:Rrs-1
        eval(['dgR',num2str(kg),'=-str2num(get(handles.r',num2str(jg),'_dg,''String''));'])
        eval(['mr',num2str(kg),'=str2num(get(handles.r',num2str(jg),'_in,''String''));']);
        if handles.g_unit==3
            eval(['Mr',num2str(kg),'=str2num(get(handles.r',num2str(jg),'_mass,''String''));']);
            eqR=eval(['[eqR,Mr',num2str(kg),'^mr',num2str(kg),'];']);
        end
        rg_store=eval(['[rg_store,mr',num2str(kg),'*dgR',num2str(kg),']']);
        mr_store=eval(['[mr_store,mr',num2str(kg),'];']);
        Symbr=eval(['{Symbr,''S',num2str(Sdn),'''};']);
        Sdn=Sdn+1;
        jg=jg+1;
    end
    handles.nP=Sdn+1;
else
    handles.nP=Sdn;
end

dg1=dg0+sum(rg_store);

%For Products

switch handles.growth
%     case 'Hoh'
%         Sdn=Sd;
% 
%         Nn=0;jh=1; pg_store=[]; m_store=[]; eqS=[]; Symbs={''};
        
    case 'Thermodynamic'
        
        m1=str2num(get(handles.p1_in,'String'));
        if handles.g_unit==3
            M1=str2num(get(handles.p1_mass,'String'));
            eqS=M1^m1;
        end
        Symbs={get(handles.p1_text,'String')}; m_store=m1;
        try
            dgP1=m1*str2num(get(handles.p1_dg,'String')); %Positive on reactant side
        catch
            dgP1=1;
        end
        jh=2; pg_store=dgP1;
end

Prs=get(handles.product_select,'Value');

if Prs>Nn
    for kg=Nn+1:Prs
        eval(['dgP',num2str(kg),'=str2num(get(handles.p',num2str(jh),'_dg,''String''));'])
        eval(['m',num2str(kg),'=str2num(get(handles.p',num2str(jh),'_in,''String''));']);
        if handles.g_unit==3
            eval(['M',num2str(kg),'=str2num(get(handles.p',num2str(jh),'_mass,''String''));']);
            eqS=eval(['[eqS,M',num2str(kg),'^m',num2str(kg),'];']);
        end
        pg_store=eval(['[pg_store,m',num2str(kg),'*dgP',num2str(kg),']']);
        m_store=eval(['[m_store,m',num2str(kg),'];']);
        Symbs=eval(['{Symbs,''S',num2str(Sdn),'''};']);
        Sdn=Sdn+1;
        jh=jh+1;
    end
end

dg2=dg1+sum(pg_store);

%Define symbolic terms
sym_Symbr=sym(Symbr); handles.symR=sym_Symbr;
sym_Symbp=sym(Symbs); handles.symP=sym_Symbp;

%Define numerator of stoichiometric substrate relationship
[ind,ind2]=find(mr_store>1);
for ll=1:sum(ind2)
    sym_Symbr(ind2(ll))=sym_Symbr(ind2(ll))^mr_store(ind2(ll));
end

[ind3,ind4]=find(m_store>1);
for mm=1:sum(ind3)
    sym_Symbp(ind4(mm))=sym_Symbp(ind4(mm))^m_store(ind4(mm));
end

SubE=prod(sym_Symbp)/prod(sym_Symbr);

if handles.g_unit==3
    MmE=log10(prod(eqR)/prod(eqS));
else
    MmE=0;
end

%Specify substrate names
%Add equations for missing substrates

%Calculate DeltaG1
RTln=R*T*log(10);

dGp = dg2+RTln*log10(SubE)+RTln*(MmE);
dG = dGp/(R*T);

%Set add substrate equation button to enabled
set(handles.add_sn,'Enable','on')
cla(handles.axes8)
axes(handles.axes8)
digits 4
text(0,0.72,['$\Delta G$ (kJ/mol) = $',char(vpa(dGp)),'$'],'interpreter','latex','horiz','left','vert','middle','fontsize',12)
digits 32
handles.dG=char(vpa(dGp,16));
handles.dG_acc=char(dGp);
switch handles.growth
    case 'Hoh'
        handles.dG=char(vpa(dG,16));
        handles.dG_acc=char(dG);
end

guidata(hObject,handles)

function dG0_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function dG0_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function r2_mass_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r2_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r3_mass_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r3_mass_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in react_sel.
function react_sel_Callback(hObject, eventdata, handles)
rsel_num=get(handles.react_sel,'Value');
csel_num=get(handles.product_select,'Value');
cla(handles.eq_disp)
if rsel_num==1
    rsel_num=2;
end

if csel_num==1;
    csel_num=2;
end
rsel_numa=rsel_num;
switch handles.motif
    case 'Cooperation'
        rsym={'S1'};
        psym={'S2'};
        num_start=3;
        rsel_num=rsel_num-2;
    case 'Threespecies'
        rsym={'S2'};
        psym={'S3'};
        
        rsel_num=rsel_num-2;
        if csel_num==1
            csel_num=2;
        end
        num_start=4;
end

for k=2:rsel_numa-1
    eval(['set(handles.r',num2str(k),'_text,''Enable'',''on'',''String'',[''S'',num2str(num_start+k-2)])'])
    eval(['set(handles.r',num2str(k),'_gas,''Enable'',''on'')'])
    if handles.g_unit<3
        eval(['set(handles.r',num2str(k),'_mass,''Enable'',''off'')'])
    else
        eval(['set(handles.r',num2str(k),'_mass,''Enable'',''on'')'])
    end
    eval(['set(handles.r',num2str(k),'_in,''Enable'',''on'')'])
end

if rsel_numa<4
    for k=rsel_numa:3
        eval(['set(handles.r',num2str(k),'_text,''Enable'',''off'')'])
        eval(['set(handles.r',num2str(k),'_gas,''Enable'',''off'')'])
        eval(['set(handles.r',num2str(k),'_mass,''Enable'',''off'')'])
        eval(['set(handles.r',num2str(k),'_in,''Enable'',''off'')'])
    end
end


rnsym={};
if rsel_num>0
    for k=num_start:num_start+rsel_num-1
        rsym=[rsym,{['+S',num2str(k)]}];
        rnsym=[rnsym,{['S',num2str(k)]}];
    end
end

pnsym={};
if csel_num>2
    z_s=k+1;
    for kk=z_s:z_s+csel_num-2
        psym=[psym,{['+S',num2str(kk)]}];
        pnsym=[pnsym,{['S',num2str(kk)]}];
    end
end


handles.react_sy=rnsym;
handles.prod_sy=pnsym;

    switch handles.motif
        case 'Cooperation'
            
            str_react=['$$',cell2mat(rsym),'\mathop{\longrightarrow}^{X1} ',cell2mat(psym),'$$'];
            axes(handles.eq_disp)
            text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
 
        case 'Threespecies'
            
            str_react=['$$',cell2mat(rsym),'\mathop{\longrightarrow}^{X2} ',cell2mat(psym),'$$'];
            axes(handles.eq_disp)
            text(0.25,0.5,str_react,'interpreter','latex','horiz','left','vert','middle','fontsize',18)
            
    end
handles.rsel_num=rsel_num;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function react_sel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in p1_gas.
function p1_gas_Callback(hObject, eventdata, handles)
p1g=get(handles.p1_gas,'Value');
if p1g==1
    set(handles.p1_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif p1g==0
    set(handles.p1_cs,'Enable','off','String','.')
    set(handles.p1_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in p2_gas.
function p2_gas_Callback(hObject, eventdata, handles)
p2g=get(handles.p2_gas,'Value');
if p2g==1
    set(handles.p2_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif p2g==0
    set(handles.p2_cs,'Enable','off','String','.')
    set(handles.p2_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in p3_gas.
function p3_gas_Callback(hObject, eventdata, handles)
p3g=get(handles.p3_gas,'Value');
if p3g==1
    set(handles.p3_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif p3g==0
    set(handles.p3_cs,'Enable','off','String','.')
    set(handles.p3_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in p4_gas.
function p4_gas_Callback(hObject, eventdata, handles)
p4g=get(handles.p4_gas,'Value');
if p4g==1
    set(handles.p4_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif p4g==0
    set(handles.p4_cs,'Enable','off','String','.')
    set(handles.p4_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
 end

% --- Executes on button press in r1_gas.
function r1_gas_Callback(hObject, eventdata, handles)
r1g=get(handles.r1_gas,'Value');
if r1g==1
    set(handles.r1_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif r1g==0
    set(handles.r1_cs,'Enable','off','String','.')
    set(handles.r1_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in r2_gas.
function r2_gas_Callback(hObject, eventdata, handles)
r2g=get(handles.r2_gas,'Value');
if r2g==1
    set(handles.r2_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif r2g==0
    set(handles.r2_cs,'Enable','off','String','.')
    set(handles.r2_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in r3_gas.
function r3_gas_Callback(hObject, eventdata, handles)
r3g=get(handles.r3_gas,'Value');
if r3g==1
    set(handles.r3_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif r3g==0
    set(handles.r3_cs,'Enable','off','String','.')
    set(handles.r3_ss,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

% --- Executes on button press in p5_gas.
function p5_gas_Callback(hObject, eventdata, handles)
p5g=get(handles.p5_gas,'Value');
if p5g==1
    set(handles.p5_ss,'Enable','on','String',handles.systems)
    set(handles.sel_comp_text,'Enable','on')
    set(handles.sel_sys_text,'Enable','on')
elseif p5g==0
    set(handles.p5_ss,'Enable','off','String','.')
    set(handles.p5_cs,'Enable','off','String','.')
    set(handles.sel_comp_text,'Enable','off')
    set(handles.sel_sys_text,'Enable','off')
end

function r3_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r2_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r2_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function r1_in_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function r1_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r1_cs.
function r1_cs_Callback(hObject, eventdata, handles)
r1v=get(handles.r1_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata1.data(r1v,:);
cmpd=get(handles.r1_cs,'String');
cmpdn=cmpd(get(handles.r1_cs,'Value'));
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_r1=polyval(p,temp);

plot(temp,dG_r1,'r.','markersize',20)
set(handles.r1_dg,'String',num2str(dG_r1),'Visible','on','Enable','on')

axes(handles.axes5)

% --- Executes during object creation, after setting all properties.
function r1_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r2_cs.
function r2_cs_Callback(hObject, eventdata, handles)
r2v=get(handles.r2_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata2.data(r2v,:);
cmpd=get(handles.r2_cs,'String');
cmpdn=cmpd(get(handles.r2_cs,'Value'));
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_r2=polyval(p,temp);

plot(temp,dG_r2,'r.','markersize',20)
set(handles.r2_dg,'String',num2str(dG_r2),'Visible','on','Enable','on')
axes(handles.axes5)


% --- Executes during object creation, after setting all properties.
function r2_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r3_cs.
function r3_cs_Callback(hObject, eventdata, handles)
r3v=get(handles.r3_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata3.data(r3v,:);
cmpd=get(handles.r3_cs,'String');
cmpdn=cmpd(get(handles.r3_cs,'Value'));
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_r3=polyval(p,temp);

plot(temp,dG_r3,'r.','markersize',20)
set(handles.r3_dg,'String',num2str(dG_r3),'Visible','on','Enable','on')
axes(handles.axes5)

% --- Executes during object creation, after setting all properties.
function r3_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rr_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function rr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p1_cs.
function p1_cs_Callback(hObject, eventdata, handles)
p1v=get(handles.p1_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata4.data(p1v,:);
cmpd=get(handles.p1_cs,'String');
cmpdn=cmpd(get(handles.p1_cs,'Value'));
set(handles.cmp1,'String',cmpdn)
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_p1=polyval(p,temp);

plot(temp,dG_p1,'r.','markersize',20)
set(handles.p1_dg,'String',num2str(dG_p1),'Visible','on','Enable','on')
axes(handles.axes5)


function p1_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p2_cs.
function p2_cs_Callback(hObject, eventdata, handles)
p2v=get(handles.p2_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata5.data(p2v,:);
cmpd=get(handles.p2_cs,'String');
cmpdn=cmpd(get(handles.p2_cs,'Value'));
set(handles.cmp2,'String',cmpdn)
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_p2=polyval(p,temp);

plot(temp,dG_p2,'r.','markersize',20)
set(handles.p2_dg,'String',num2str(dG_p2),'Visible','on','Enable','on')
axes(handles.axes5)

% --- Executes during object creation, after setting all properties.
function p2_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p3_cs.
function p3_cs_Callback(hObject, eventdata, handles)
p3v=get(handles.p3_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata6.data(p3v,:);
cmpd=get(handles.p3_cs,'String');
cmpdn=cmpd(get(handles.p3_cs,'Value'));
set(handles.cmp3,'String',cmpdn)
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_p3=polyval(p,temp);

plot(temp,dG_p3,'r.','markersize',20)
set(handles.p3_dg,'String',num2str(dG_p3),'Visible','on','Enable','on')
axes(handles.axes5)

% --- Executes during object creation, after setting all properties.
function p3_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p4_cs.
function p4_cs_Callback(hObject, eventdata, handles)
p4v=get(handles.p4_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata7.data(p4v,:);
cmpd=get(handles.p4_cs,'String');
cmpdn=cmpd(get(handles.p4_cs,'Value'));
set(handles.cmp4,'String',cmpdn)
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_p4=polyval(p,temp);

plot(temp,dG_p4,'r.','markersize',20)

set(handles.p4_dg,'String',num2str(dG_p4),'Visible','on','Enable','on')
axes(handles.axes5)
% --- Executes during object creation, after setting all properties.
function p4_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p5_cs.
function p5_cs_Callback(hObject, eventdata, handles)
p5v=get(handles.p5_cs,'Value');
temp=str2num(get(handles.temp_in,'String'));
compound = handles.Gdata8.data(p5v,:);
cmpd=get(handles.p5_cs,'String');
cmpdn=cmpd(get(handles.p5_cs,'Value'));
set(handles.cmp5,'String',cmpdn)
axes(handles.dG_ax)
hold off
plot(handles.temperatures,compound,'linewidth',2)
xlabel('Temperature (K)')
ylabel('\Delta G^0 (kJ/mol)')
legend(cmpdn)
hold on
%Calculate fit
p = polyfit(handles.temperatures,compound,2);
dG_p5=polyval(p,temp);

plot(temp,dG_p5,'r.','markersize',20)
set(handles.p5_dg,'String',num2str(dG_p5),'Visible','on','Enable','on')
axes(handles.axes5)

% --- Executes during object creation, after setting all properties.
function p5_cs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r1_ss.
function r1_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.r1_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata1=load('system1.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata1=load('system2.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata1=load('system3.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata1=load('system4.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata1=load('system5.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata1=load('system6.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata1=load('system7.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata1=load('system8.mat','-mat');
        set(handles.r1_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function r1_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r2_ss.
function r2_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.r2_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata2=load('system1.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata2=load('system2.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata2=load('system3.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata2=load('system4.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata2=load('system5.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata2=load('system6.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata2=load('system7.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata2=load('system8.mat','-mat');
        set(handles.r2_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
function r2_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in r3_ss.
function r3_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.r3_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata3=load('system1.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata3=load('system2.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata3=load('system3.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata3=load('system4.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata3=load('system5.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata3=load('system6.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata3=load('system7.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata3=load('system8.mat','-mat');
        set(handles.r3_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function r3_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p1_ss.
function p1_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.p1_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata4=load('system1.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata4=load('system2.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata4=load('system3.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata4=load('system4.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata4=load('system5.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata4=load('system6.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata4=load('system7.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata4=load('system8.mat','-mat');
        set(handles.p1_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function p1_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p2_ss.
function p2_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.p2_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata5=load('system1.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata5=load('system2.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata5=load('system3.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata5=load('system4.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata5=load('system5.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata5=load('system6.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata5=load('system7.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata5=load('system8.mat','-mat');
        set(handles.p2_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function p2_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p3_ss.
function p3_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.p3_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata6=load('system1.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata6=load('system2.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata6=load('system3.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata6=load('system4.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata6=load('system5.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata6=load('system6.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata6=load('system7.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata6=load('system8.mat','-mat');
        set(handles.p3_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function p3_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p4_ss.
function p4_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.p4_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata7=load('system1.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata7=load('system2.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata7=load('system3.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata7=load('system4.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata7=load('system5.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata7=load('system6.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata7=load('system7.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata7=load('system8.mat','-mat');
        set(handles.p4_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function p4_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in p5_ss.
function p5_ss_Callback(hObject, eventdata, handles)
ss_v=get(handles.p5_ss,'Value');

sel_sys=handles.systems{ss_v};

switch sel_sys
    case 'H-O'
        cd('thermodyamic_data')
        handles.Gdata8=load('system1.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds1)
        cd ..
    case 'H-O-N'
        cd('thermodyamic_data')
        handles.Gdata8=load('system2.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds2)
        cd ..
    case 'H-O-S'
        cd('thermodyamic_data')
        handles.Gdata8=load('system3.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds3)
        cd ..
    case 'H-O-N-S-C_inorganic'
        cd('thermodyamic_data')
        handles.Gdata8=load('system4.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds4)
        cd ..
    case 'H-O-C_organic'
        cd('thermodyamic_data')
        handles.Gdata8=load('system5.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds5)
        cd ..
    case 'H-O-S-C_metals/minerals'
        cd('thermodyamic_data')
        handles.Gdata8=load('system6.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds6)
        cd ..
    case 'H-O-P'
        cd('thermodyamic_data')
        handles.Gdata8=load('system7.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds7)
        cd ..
    case 'H-O-Halogens'
        cd('thermodyamic_data')
        handles.Gdata8=load('system8.mat','-mat');
        set(handles.p5_cs,'Enable','on','String',handles.compounds8)
        cd ..
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function p5_ss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function r1_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function r2_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function r3_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function p1_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function p2_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function p3_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function p4_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function p5_dg_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p1_gamma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p1_gamma_CreateFcn(hObject, eventdata, handles)


function p2_gamma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p2_gamma_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p3_gamma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p3_gamma_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p4_gamma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p4_gamma_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function p5_gamma_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function p5_gamma_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in add_sn.
function add_sn_Callback(hObject, eventdata, handles)
%Define equation for additional reactants
handles.symR;

%Define equations for additional substrates (Products)
syms D
noP=length(handles.symP); gamma1=[];

switch handles.growth
    
    case 'Thermodynamic'

        P=[]; Pt=[];
        if noP>1
            for k=1:noP-1
                gamma1{k}=get(eval(['handles.p',num2str(k+length(handles.gammas1)),'_gamma']),'String');
                eq{k}=-D*handles.symP(k+1)+gamma1{k}*handles.funcN;
                eqtx{k}=['$\frac{dS_',num2str(handles.nP+(k-1)),'}{dt} = -DS_',num2str(handles.nP+(k-1)),handles.funcTa,num2str(k+2+(length(handles.gammas1)-1)),handles.funcTb,'$'];
                eval(['eq',num2str(handles.eqno+1),'=eq{k};'])
                P=[P,eq(k)];
                Pt=[Pt;eqtx{k}];
            end
        else
            gamma1=get(handles.p1_gamma,'String');
            eq=-D*handles.symP+gamma1*handles.funcN;
            eqtx=['$\frac{dS_',num2str(handles.nP),'}{dt} = -DS_',num2str(handles.nP),handles.funcTa,num2str(handles.rte+1),handles.funcTb,'$'];
            eval(['eq',num2str(handles.eqno+1),'=eq;'])
            P=[P,eq];
            Pt=[Pt;eqtx];
        end
        Ptn=Pt;
        
%     case 'Hoh'
%         if handles.xnum>1
%             
%             cmp_names2=[];
%             for k=1:length(handles.prod_sy)
%                 compN=eval(['get(handles.cmp',num2str(k),',''String'');']);
%                 if isempty(compN)
%                     compN=['NewN',num2str(k)];
%                 end
%                 cmp_names2=strvcat(cmp_names2,compN);
%             end
%             [inter_names,inter_indx]=intersect(cmp_names2,handles.old_cmps); %Find location of matching products in both reactions (in new)
%             [inter_names_o,inter_indx_o]=intersect(handles.old_cmps,cmp_names2); %Find location of matching products in both reactions (in old)
%             
%             [inter_names2,inter_indx2]=setdiff(cmp_names2,handles.old_cmps); %Find location of unique products in both reactions (in new)
%             [inter_names2_o,inter_indx2_o]=setdiff(handles.old_cmps,cmp_names2); %Find location of unique products in both reactions (in old)
%             eq=handles.old_eqs{1}; eqtx=handles.old_eqtx;
%             
%             P=[]; Pt=[]; Prod=[]; eq1={}; eq2={}; eqtx1=[];
%             %Assign matching compounds to same equation
%             if ~isempty(inter_indx)
%                 for kk=1:length(inter_indx)
%                     jk=inter_indx(kk);
%                     jko=inter_indx_o(kk);
%                     gamma1{jko}=get(eval(['handles.p',num2str(jk),'_gamma']),'String');
%                     eq1{jko}=eq{jko}+gamma1{jko}*handles.funcN;
%                     eqtx1{jko}=[eqtx(jko,1:end-1),handles.funcTa,num2str(handles.gamma_count(end)+jk),handles.funcTb,'$'];
%                     eval(['eq',num2str(handles.eqno+kk),'=eq;'])
%                     P=[P,eq1{jko}];
%                     Pt=strvcat(Pt,eqtx1{jko});
%                     Prod=[Prod,{eval(['get(handles.p',num2str(jko),'_text,''String'')'])}];
%                 end
%             end
%             
%             eq1a=eq(inter_indx2_o); %Add old unique equations
%             P=[P,eq1a];
% 
%          %Assign unique compounds to new equation
%          if ~isempty(inter_indx2)
%          un_indx=1+length(handles.old_cmps)+length(cmp_names2)-length(inter_indx);
%          gam_indx=handles.gamma_count(end)+kk;
%            for kj=1:length(inter_indx2)
%                 jj=inter_indx2(kj);
%                 jjo=inter_indx2_o(kj);
%                 gamma1{jjo}=get(eval(['handles.p',num2str(jj),'_gamma']),'String');
%                 eq2{jjo}=-D*sym(['S',num2str(un_indx)])+gamma1{jjo}*handles.funcN;
%                 eqtx=['$\frac{dS_',num2str(un_indx),'}{dt} = -DS_',num2str(un_indx),handles.funcTa,num2str(handles.gamma_count(end)+jj),handles.funcTb,'$'];
%                 eval(['eq',num2str(handles.eqno+kj),'=eq;'])
%                 P=[P,eq2];
%                 Pt=strvcat(Pt,eqtx);
%                 Prod=[Prod,{['S',num2str(un_indx)]}];
% 
%            end
%            
%            %Add previous equations
%            Ptn=strvcat(handles.old_eqtx(inter_indx2_o,:),Pt);
%            
%            %Remove empty cells
%            try
%                P = P(~cellfun('isempty',P));
%            end
%            %Set reaction equation equal to correct substrates and change
%            %variable names and dG
%            var_diff1=setdiff(handles.prod_sy,Prod);
%            var_diff2=setdiff(Prod,handles.prod_sy);
%            
%            handles.prod_sy=Prod;
%            handles.dG=subs(handles.dG,var_diff1,var_diff2);
%            handles.dG_acc=subs(handles.dG_acc,var_diff1,var_diff2);
%            cla(handles.axes8)
%            axes(handles.axes8)
%            digits 4
%            text(0,0.72,['$\Delta G$ (kJ/mol) = $',char(vpa(handles.dG_acc)),'$'],'interpreter','latex','horiz','left','vert','middle','fontsize',12)
%            digits 32
%          else
%              Ptn=Pt;
%          end
%             
%         else
%             
%             P=[]; Pt=[]; handles.gamma_count=[];
%             if noP>1
%                 for k=1:noP
%                     gamma1{k}=get(eval(['handles.p',num2str(k+length(handles.gammas1)),'_gamma']),'String');
%                     eq{k}=-D*handles.symP(k)+gamma1{k}*handles.funcN;
%                     eqtx{k}=['$\frac{dS_',num2str(handles.nP+(k-2)),'}{dt} = -DS_',num2str(handles.nP+(k-2)),handles.funcTa,num2str(k+1+(length(handles.gammas1)-1)),handles.funcTb,'$'];
%                     handles.gamma_count=[handles.gamma_count,k+1+(length(handles.gammas1)-1)];
%                     eval(['eq',num2str(handles.eqno+1),'=eq{k};'])
%                     P=[P,eq(k)];
%                     Pt=[Pt;eqtx{k}];
%                 end
%             else
%                 gamma1=get(handles.p1_gamma,'String');
%                 eq=-D*handles.symP+gamma1*handles.funcN;
%                 eqtx=['$\frac{dS_',num2str(handles.nP),'}{dt} = -DS_',num2str(handles.nP),handles.funcTa,num2str(handles.rte+1),handles.funcTb,'$'];
%                 eval(['eq',num2str(handles.eqno+1),'=eq;'])
%                 P=[P,eq];
%                 Pt=[Pt;eqtx];
%             end
%             Ptn=Pt;
%         end
end

for k=1:noP
    gamma{k}=get(eval(['handles.p',num2str(k),'_gamma']),'String');
end
cla(handles.axes7)
axes(handles.axes7)
text(0.3,0.4,Pt,'interpreter','latex','FontSize',14)
set(handles.finish_but,'Enable','on')
handles.thermo_eqs=P;
handles.thermo_text=Ptn;
handles.gamma_out=gamma;
guidata(hObject,handles)

% --- Executes on button press in finish_but.
function finish_but_Callback(hObject, eventdata, handles)
handles.output1=handles.thermo_eqs;
handles.output2=handles.thermo_text;
handles.output3=handles.gamma_out;
handles.output4=[handles.react_sy,handles.prod_sy];
%handles.kr=str2num(get(handles.rr,'String'));
% if ~isnumeric(handles.kr) || handles.kr>1
%     handles.kr=1;
% end

if handles.xnum==1
    handles.cmp_names=[];
    for k=1:length(handles.prod_sy)
        
        comp=eval(['get(handles.cmp',num2str(k),',''String'');']);
        if isempty(comp)
            comp=['New',num2str(k)];
        end
        handles.cmp_names=strvcat(handles.cmp_names,comp);
        
    end
end
   
guidata(hObject,handles)
uiresume(handles.figure1);

function cmp1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cmp1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cmp2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cmp2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cmp3_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cmp3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cmp4_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cmp4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cmp5_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function cmp5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in unit_sel.
function unit_sel_Callback(hObject, eventdata, handles)
handles.g_unit=get(handles.unit_sel,'Value');
if handles.g_unit < 3
    set(handles.mass_text,'Enable','off') %Molar units
else
    set(handles.mass_text,'Enable','on') %COD units
    set(handles.r1_mass,'Enable','on')
    set(handles.p1_mass,'Enable','on')
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function unit_sel_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
