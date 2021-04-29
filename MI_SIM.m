function varargout = MI_SIM(varargin)
% MI_SIM MATLAB code for MI_SIM.fig
%      MI_SIM, by itself, creates a new MI_SIM or raises the existing
%      singleton*.
%
%      H = MI_SIM returns the handle to a new MI_SIM or the handle to
%      the existing singleton*.
%
%      MI_SIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MI_SIM.M with the given input arguments.
%
%      MI_SIM('Property','Value',...) creates a new MI_SIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MI_SIM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MI_SIM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MI_SIM

% Last Modified by GUIDE v2.5 24-Jan-2017 11:38:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MI_SIM_OpeningFcn, ...
    'gui_OutputFcn',  @MI_SIM_OutputFcn, ...
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


% --- Executes just before MI_SIM is made visible.
function MI_SIM_OpeningFcn(hObject, eventdata, handles, varargin)
%pwd=MISIM_nu2016;

%reset the plots
axes(handles.solutionplot);
cla reset
axes(handles.trajectoryplot);
cla reset

xlabel(handles.solutionplot,'Time (days)')
ylabel(handles.solutionplot,'Concentration (kgCOD m^{-3})')

set(handles.twodimplot,'Checked','on'); set(handles.threedimplot,'Checked','off')
set(handles.timeplot,'Checked','on'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','off')
set(handles.sol_disp,'String','Display:')

set(handles.overlay,'enable','off')

%reset all of the values

set(handles.solver,'Value',3);
set(handles.abstol,'String','1e-8');
set(handles.reltol,'String','1e-8');
set(handles.step_size,'String','0.01');
set(handles.prec_num,'String','32');
set(handles.error_val,'String','1e-6');
set(handles.pert_p,'String','0.0001');
set(handles.func_prog,'String','');
set(handles.lsanaly,'Value',1,'ForegroundColor',[0.078, 0.169, 0.549])
set(handles.routhcrit,'Value',0,'ForegroundColor','r')
set(handles.jacobian_but,'Value',0,'ForegroundColor','r')
set(handles.uipanel6,'Title','Plot of trajectory from initial conditions')
set(handles.spmatrix,'Visible','off');
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off');
set(handles.normeig,'Visible','off','Value',0)
set(handles.infotable,'Visible','off')
set(handles.stabtable,'Visible','off')
set(handles.gamma3,'Visible','off','enable','off'); set(handles.gam_text3,'Visible','off','enable','off')
set(handles.gamma4,'Visible','off','enable','off'); set(handles.gam_text4,'Visible','off','enable','off')
set(handles.gamma5,'Visible','off','enable','off'); set(handles.gam_text5,'Visible','off','enable','off')
set(handles.gamma6,'Visible','off','enable','off'); set(handles.gam_text6,'Visible','off','enable','off')

set(handles.furth_plots,'Visible','off')
handles.simtype='single_p';
handles.growthmodel='Monod';
handles.gui_front=0;
%For thermodynamics
handles.gammas=[];
handles.thermeqs=[];
handles.dG=[];
handles.Temperature=[];

handles.steqs=[]; %Extra variables if thermodynamics used

dir_rep=dir('Reports');

if length(dir_rep)<5
    set(handles.sendreport,'Enable','off')
    set(handles.report_menu,'Enable','off')
    set(handles.reportpanel,'Visible','off'); set(handles.email_add,'Enable','off'); set(handles.server_add,'Enable','off')
    set(handles.email_txt,'Enable','off'); set(handles.serv_txt,'Enable','off');
else
    set(handles.report_menu,'Enable','on')
    set(handles.sendreport,'Enable','on')
    set(handles.reportpanel,'Visible','on');set(handles.email_add,'Enable','on'); set(handles.server_add,'Enable','on')
    set(handles.email_txt,'Enable','on'); set(handles.serv_txt,'Enable','on');
end
handles.yout_phase=[];
handles.clcyc=1;

%List of Existing Growth Models
Exist_Models=[{'Commensalism'};{'Competition'};{'Predation'};{'No_interaction'};{'Cooperation'};{'Amensalism'};{'Threespecies'}];
Exist_Gmodels=[{'Monod'};{'Contois'};{'Tessier'};{'Moser'};{'Haldane'};{'Andrews'};{'Thermodynamic'}];
if ~isempty(varargin)
    %Check valid model
    if numel(varargin)==2
        fndmd=strncmpi(varargin(1),Exist_Models,4); %Check first four letters of model for match
        fndmg=strncmpi(varargin(2),Exist_Gmodels,3); %Check first three letters of model for match
    else
        try
            fndmd=strncmpi(varargin,Exist_Models,4);
            fndmg=0;
        catch
            fndmg=strncmpi(varargin(2),Exist_Gmodels,3);
            fndmd=0;
        end
    end
    mdindex=find(fndmd); %Get index
    mgindex=find(fndmg);
    if (sum(fndmd)>0)
        handles.motif_name=Exist_Models{mdindex};
    else
        handles.motif_name='Cooperation';
    end
    
    if (sum(fndmg)>0)
        handles.growthmodel=Exist_Gmodels{mgindex};
        set(handles.growthmenu,'Value',mgindex);
    else
        handles.growthmodel='Monod';
        set(handles.growthmenu,'Value',1);
    end
else
    handles.motif_name='Syntrophy'; handles.motif='syn';
end
try
    handles=model_sel(handles);
catch
    return
end

% Choose default command line output for MI_SIM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MI_SIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MI_SIM_OutputFcn(hObject, eventdata, handles)
try
    varargout{1} = handles.output;
catch
    return
end

function x1_init_Callback(hObject, eventdata, handles)
X1_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function x1_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x2_init_Callback(hObject, eventdata, handles)
X2_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function x2_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s1_init_Callback(hObject, eventdata, handles)
S1_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s1_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function s2_init_Callback(hObject, eventdata, handles)
S2_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s2_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function s3_init_Callback(hObject, eventdata, handles)
S3_init=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s3_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x3_init_Callback(hObject, eventdata, handles)
X3_init=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function x3_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s4_init_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function s4_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s5_init_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function s5_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s6_init_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function s6_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in growthmenu.
function growthmenu_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.growthmodel=contents{get(hObject,'Value')};

handles=model_sel(handles);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function growthmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km1_in_Callback(hObject, eventdata, handles)

km1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function km1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y1_in_Callback(hObject, eventdata, handles)
Y1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function y1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec1_in_Callback(hObject, eventdata, handles)
kdec1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function kdec1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km2_in_Callback(hObject, eventdata, handles)
km2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function km2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y2_in_Callback(hObject, eventdata, handles)
Y2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function y2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec2_in_Callback(hObject, eventdata, handles)
kdec2=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function kdec2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ki2_in_Callback(hObject, eventdata, handles)
KI2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function ki2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks1_in_Callback(hObject, eventdata, handles)
Ks1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks2_in_Callback(hObject, eventdata, handles)
Ks2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function ks2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks3_in_Callback(hObject, eventdata, handles)
Ks3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n1_in_Callback(hObject, eventdata, handles)
n1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function n1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n2_in_Callback(hObject, eventdata, handles)
n2=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function n2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km3_in_Callback(hObject, eventdata, handles)
km3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function km3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y3_in_Callback(hObject, eventdata, handles)
Y3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function y3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec3_in_Callback(hObject, eventdata, handles)
kdec3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function kdec3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function n3_in_Callback(hObject, eventdata, handles)
n3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function n3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function d_in_Callback(hObject, eventdata, handles)
D=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function d_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s1in_in_Callback(hObject, eventdata, handles)
S1in=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s1in_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s2in_in_Callback(hObject, eventdata, handles)
S2in=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s2in_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s3in_in_Callback(hObject, eventdata, handles)
S3in=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s3in_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function time_in_Callback(hObject, eventdata, handles)
time1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function time_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma0_Callback(hObject, eventdata, handles)
gamma0=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ks32_in_Callback(hObject, eventdata, handles)
ks3c=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks32_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma1_Callback(hObject, eventdata, handles)
gamma1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma2_Callback(hObject, eventdata, handles)
gamma2=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma3_Callback(hObject, eventdata, handles)
gamma3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma4_Callback(hObject, eventdata, handles)
gamma4=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma4_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma5_Callback(hObject, eventdata, handles)
gamma5=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma6_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function gamma6_CreateFcn(hObject, eventdata, handles)
gamma6=str2double(get(hObject,'String'));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in solver.
function solver_Callback(hObject, eventdata, handles)
solvstr=get(hObject,'String');
solvnum=get(hObject,'Value');
solver=solvstr(solvnum,:);
jacb=get(handles.jacobian_but,'Value');
switch strtrim(char(solver))
    case 'ode23s'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode15s'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode23t'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode23'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
    case 'ode45'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
    case 'ode113'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
end

% --- Executes during object creation, after setting all properties.
function solver_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function abstol_Callback(hObject, eventdata, handles)
abstol=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function abstol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reltol_Callback(hObject, eventdata, handles)
reltol=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function reltol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function error_val_Callback(hObject, eventdata, handles)
error=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function error_val_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function step_size_Callback(hObject, eventdata, handles)
stepsize=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function step_size_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function prec_num_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function prec_num_CreateFcn(hObject, eventdata, handles)
precision=str2double(get(hObject,'string'));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_button2.
function plot_button2_Callback(hObject, eventdata, handles)
which_plot='traj';
plot_results;

% --- Executes on button press in s1_check.
function s1_check_Callback(hObject, eventdata, handles)
s1_check=get(hObject,'Value');


% --- Executes on button press in x1_check.
function x1_check_Callback(hObject, eventdata, handles)
x1_check=get(hObject,'Value');


% --- Executes on button press in s2_check.
function s2_check_Callback(hObject, eventdata, handles)
s2_check=get(hObject,'Value');


% --- Executes on button press in x2_check.
function x2_check_Callback(hObject, eventdata, handles)
x2_check=get(hObject,'Value');

% --- Executes on button press in s3_check.
function s3_check_Callback(hObject, eventdata, handles)
s3_check=get(hObject,'Value');

% --- Executes on button press in x3_check.
function x3_check_Callback(hObject, eventdata, handles)
x3_check=get(hObject,'Value');

% --- Executes on button press in s4_check.
function s4_check_Callback(hObject, eventdata, handles)
s4_check=get(hObject,'Value');

% --- Executes on button press in s5_check.
function s5_check_Callback(hObject, eventdata, handles)
s5_check=get(hObject,'Value');

% --- Executes on button press in s6_check.
function s6_check_Callback(hObject, eventdata, handles)
s6_check=get(hObject,'Value');

% --- Executes on button press in overlay.
function overlay_Callback(hObject, eventdata, handles)
multi=get(hObject,'Value');

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

% Get user input
growth=handles.growthmodel;
motif=handles.motif;

km1 = str2double(get(handles.km1_in,'String'));
if km1<0 || isnumeric(km1)==0
    set(handles.km1_in,'String',13);
    km1 = str2double(get(handles.km1_in,'String'));
    msgbox('The value entered for km1 was not valid so the default value was reset','Error','error')
else
end
Y1 = str2double(get(handles.y1_in,'String'));
if Y1<0 || isnumeric(Y1)==0
    set(handles.y1_in,'String',0.04);
    Y1 = str2double(get(handles.y1_in,'String'));
    msgbox('The value entered for Y1 was not valid so the default value was reset','Error','error')
else
end
kdec1 = str2double(get(handles.kdec1_in,'String'));
if kdec1<0 || kdec1>1 || isnumeric(kdec1)==0
    set(handles.kdec1_in,'String',0.02);
    kdec1 = str2double(get(handles.kdec1_in,'String'));
    msgbox('The value entered for (0<)kdec1(<1) was not valid so the default value was reset','Error','error')
else
end
km2 = str2double(get(handles.km2_in,'String'));
if km2<0 || isnumeric(km2)==0
    set(handles.km2_in,'String',35);
    km2 = str2double(get(handles.km2_in,'String'));
    msgbox('The value entered for km2 was not valid so the default value was reset','Error','error')
else
end
Y2 = str2double(get(handles.y2_in,'String'));
if Y2<0 || isnumeric(Y2)==0
    set(handles.y2_in,'String',0.06);
    Y2 = str2double(get(handles.y2_in,'String'));
    msgbox('The value entered for Y2 was not valid so the default value was reset','Error','error')
else
end
kdec2 = str2double(get(handles.kdec2_in,'String'));
if kdec2<0 || kdec2>1 || isnumeric(kdec2)==0
    set(handles.kdec2_in,'String',0.02);
    kdec2 = str2double(get(handles.kdec2_in,'String'));
    msgbox('The value entered for (0<)kdec2(<1) was not valid so the default value was reset','Error','error')
else
end

km3 = str2double(get(handles.km3_in,'String'));
if km3<0 || isnumeric(km3)==0
    set(handles.km3_in,'String',21);
    km3 = str2double(get(handles.km3_in,'String'));
    msgbox('The value entered for km3 was not valid so the default value was reset','Error','error')
else
end
Y3 = str2double(get(handles.y3_in,'String'));
if Y3<0 || isnumeric(Y3)==0
    set(handles.y3_in,'String',0.04);
    Y3 = str2double(get(handles.y3_in,'String'));
    msgbox('The value entered for Y3 was not valid so the default value was reset','Error','error')
else
end
kdec3 = str2double(get(handles.kdec3_in,'String'));
if kdec3<0 || kdec3>1 || isnumeric(kdec3)==0
    set(handles.kdec3_in,'String',0.02);
    kdec3 = str2double(get(handles.kdec3_in,'String'));
    msgbox('The value entered for (0<)kdec3(<1) was not valid so the default value was reset','Error','error')
else
end

S1in = str2double(get(handles.s1in_in,'String'));
if S1in<0 || isnumeric(S1in)==0
    set(handles.s1in_in,'String',5);
    S1in = str2double(get(handles.s1in_in,'String'));
    msgbox('The value entered for S1in was not valid so the default value was reset','Error','error')
else
end
S2in = str2double(get(handles.s2in_in,'String'));
if S2in<0 || isnumeric(S2in)==0
    set(handles.s2in_in,'String',5);
    S2in = str2double(get(handles.s2in_in,'String'));
    msgbox('The value entered for S2in was not valid so the default value was reset','Error','error')
else
end

S3in = str2double(get(handles.s3in_in,'String'));
if S3in<0 || isnumeric(S3in)==0
    set(handles.s3in_in,'String',1);
    S3in = str2double(get(handles.s3in_in,'String'));
    msgbox('The value entered for S3in was not valid so the default value was reset','Error','error')
else
end
KI2 = str2double(get(handles.ki2_in,'String'));
if KI2<0 || isnumeric(KI2)==0
    set(handles.ki2_in,'String',0.0000035);
    KI2 = str2double(get(handles.ki2_in,'String'));
    msgbox('The value entered for KI was not valid so the default value was reset','Error','error')
else
end
Ks1 = str2double(get(handles.ks1_in,'String'));
if Ks1<0 || isnumeric(Ks1)==0
    set(handles.ks1_in,'String',0.3);
    Ks1 = str2double(get(handles.ks1_in,'String'));
    msgbox('The value entered for KS1) was not valid so the default value was reset','Error','error')
else
end
Ks2 = str2double(get(handles.ks2_in,'String'));
if Ks2<0 || isnumeric(Ks2)==0
    set(handles.ks2_in,'String',0.000025);
    Ks2 = str2double(get(handles.ks2_in,'String'));
    msgbox('The value entered for KS2 was not valid so the default value was reset','Error','error')
else
end

Ks3 = str2double(get(handles.ks3_in,'String'));
if Ks3<0 || isnumeric(Ks3)==0
    set(handles.ks3_in,'String',0.001);
    Ks3 = str2double(get(handles.ks3_in,'String'));
    msgbox('The value entered for KS3 was not valid so the default value was reset','Error','error')
else
end

Ks3c = str2double(get(handles.ks32_in,'String'));
if Ks3c<0 || isnumeric(Ks3c)==0
    set(handles.ks32_in,'String',1e-6);
    Ks3c = str2double(get(handles.ks32_in,'String'));
    msgbox('The value entered for KS3c was not valid so the default value was reset','Error','error')
else
end

n1 = str2double(get(handles.n1_in,'String'));
if n1<0 || isnumeric(n1)==0
    set(handles.n1_in,'String',3);
    n1 = str2double(get(handles.n1_in,'String'));
    msgbox('The value entered for n1 was not valid so the default value was reset','Error','error')
else
end
n2 = str2double(get(handles.n2_in,'String'));
if n2<0 || isnumeric(n2)==0
    set(handles.n2_in,'String',2);
    n2 = str2double(get(handles.n2_in,'String'));
    msgbox('The value entered for n2 was not valid so the default value was reset','Error','error')
else
end
n3 = str2double(get(handles.n3_in,'String'));
if n3<0 || isnumeric(n3)==0
    set(handles.n3_in,'String',2);
    n3 = str2double(get(handles.n3_in,'String'));
    msgbox('The value entered for n3 was not valid so the default value was reset','Error','error')
else
end
gamma0 = str2double(get(handles.gamma0,'String'));
if gamma0<0 || isnumeric(gamma0)==0
    set(handles.gamma0,'String',0.43);
    gamma0 = str2double(get(handles.gamma0,'String'));
    msgbox('The value entered for \gamma_0 was not valid so the default value was reset','Error','error')
else
end
gamma1 = str2double(get(handles.gamma1,'String'));
if gamma1<0 || isnumeric(gamma1)==0
    set(handles.gamma1,'String',0.1429);
    gamma1 = str2double(get(handles.gamma1,'String'));
    msgbox('The value entered for \gamma_1 was not valid so the default value was reset','Error','error')
else
end
gamma2 = str2double(get(handles.gamma2,'String'));
if gamma2<0 || isnumeric(gamma2)==0
    set(handles.gamma2,'String',0.0769);
    gamma2 = str2double(get(handles.gamma2,'String'));
    msgbox('The value entered for \gamma_2 was not valid so the default value was reset','Error','error')
else
end
D = str2double(get(handles.d_in,'String'));
if D<0 || isnumeric(D)==0
    set(handles.d_in,'String',0.1);
    D = str2double(get(handles.d_in,'String'));
    msgbox('The value entered for dilution(>0) was not valid so the default value was reset','Error','error')
else
end
time1 = str2double(get(handles.time_in,'String'));
if time1<0 || isnumeric(time1)==0
    set(handles.time_in,'String',1000);
    time1 = str2double(get(handles.time_in,'String'));
    msgbox('The value entered for time(>0) was not valid so the default value was reset','Error','error')
else
end

S1_init=str2double(get(handles.s1_init,'String'));
if S1_init<0 || isnumeric(S1_init)==0
    set(handles.s1_init,'String',0.1);
    S1_init=str2double(get(handles.s1_init,'String'));
    msgbox('The value entered for S1_init was not valid so the default value was reset','Error','error')
else
end
X1_init=str2double(get(handles.x1_init,'String'));
if X1_init<0 || isnumeric(X1_init)==0
    set(handles.x1_init,'String',0.1);
    X1_init=str2double(get(handles.x1_init,'String'));
    msgbox('The value entered for X1_init was not valid so the default value was reset','Error','error')
else
end
S2_init=str2double(get(handles.s2_init,'String'));
if S2_init<0 || isnumeric(S2_init)==0
    set(handles.s2_init,'String',0.1);
    S2_init=str2double(get(handles.s2_init,'String'));
    msgbox('The value entered for S2_init was not valid so the default value was reset','Error','error')
else
end
X2_init=str2double(get(handles.x2_init,'String'));
if X2_init<0 || isnumeric(X2_init)==0
    set(handles.x2_init,'String',0.1);
    X2_init=str2double(get(handles.x2_init,'String'));
    msgbox('The value entered for X2_init was not valid so the default value was reset','Error','error')
else
end
S3_init=str2double(get(handles.s3_init,'String'));
if S3_init<0 || isnumeric(S3_init)==0
    set(handles.s3_init,'String',0.1);
    S3_init=str2double(get(handles.s3_init,'String'));
    msgbox('The value entered for S3_init was not valid so the default value was reset','Error','error')
else
end
X3_init=str2double(get(handles.x3_init,'String'));
if X3_init<0 || isnumeric(X3_init)==0
    set(handles.x3_init,'String',0.1);
    X3_init=str2double(get(handles.x3_init,'String'));
    msgbox('The value entered for X3_init was not valid so the default value was reset','Error','error')
else
end

S4_init=str2double(get(handles.s4_init,'String'));
if S4_init<0 || isnumeric(S4_init)==0
    set(handles.s4_init,'String',0.1);
    S4_init=str2double(get(handles.s4_init,'String'));
    msgbox('The value entered for S4_init was not valid so the default value was reset','Error','error')
else
end

S5_init=str2double(get(handles.s5_init,'String'));
if S5_init<0 || isnumeric(S3_init)==0
    set(handles.s5_init,'String',0.1);
    S5_init=str2double(get(handles.s5_init,'String'));
    msgbox('The value entered for S5_init was not valid so the default value was reset','Error','error')
else
end

S6_init=str2double(get(handles.s6_init,'String'));
if S6_init<0 || isnumeric(S6_init)==0
    set(handles.s6_init,'String',0.1);
    S6_init=str2double(get(handles.s6_init,'String'));
    msgbox('The value entered for S6_init was not valid so the default value was reset','Error','error')
else
end

abstol=str2double(get(handles.abstol,'String'));
if abstol<0 || isnumeric(abstol)==0
    set(handles.abstol,'String',1e-8);
    abstol=str2double(get(handles.abstol,'String'));
    msgbox('The value entered for Abs. tolerance was not valid so the default value was reset','Error','error')
else
end

reltol=str2double(get(handles.reltol,'String'));
if reltol<0 || isnumeric(reltol)==0
    set(handles.reltol,'String',1e-8);
    reltol=str2double(get(handles.reltol,'String'));
    msgbox('The value entered for Rel. tolerance was not valid so the default value was reset','Error','error')
else
end

stepsize=str2double(get(handles.step_size,'String'));
if stepsize<0 || isnumeric(stepsize)==0
    set(handles.step_size,'String',0.01);
    stepsize=str2double(get(handles.step_size,'String'));
    msgbox('The value entered for step size was not valid so the default value was reset','Error','error')
else
end

solvstr=get(handles.solver,'String');
solvnum=get(handles.solver,'Value');
if solvnum==1
    solvnum=3;
end
solver=solvstr(solvnum,:);

s1_check=get(handles.s1_check,'Value');
x1_check=get(handles.x1_check,'Value');
s2_check=get(handles.s2_check,'Value');
x2_check=get(handles.x2_check,'Value');
s3_check=get(handles.s3_check,'Value');
x3_check=get(handles.x3_check,'Value');
s4_check=get(handles.s4_check,'Value');
s5_check=get(handles.s5_check,'Value');
s6_check=get(handles.s6_check,'Value');

run_script_gui;

function fixed_points_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function fixed_points_s1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function models_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function commensalism_Callback(hObject, eventdata, handles)
handles.motif_name='Commensalism';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};
handles=model_sel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function competition_Callback(hObject, eventdata, handles)
handles.motif_name='Competition';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function predation_Callback(hObject, eventdata, handles)
handles.motif_name='Predation';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function no_interaction_Callback(hObject, eventdata, handles)
handles.motif_name='No_interaction';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);

guidata(hObject,handles)

% --------------------------------------------------------------------
function cooperation_Callback(hObject, eventdata, handles)
handles.motif_name='Cooperation';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);

guidata(hObject,handles)

% --------------------------------------------------------------------
function amensalism_Callback(hObject, eventdata, handles)
handles.motif_name='Amensalism';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);
guidata(hObject,handles)

% --------------------------------------------------------------------
function threesp_Callback(hObject, eventdata, handles)
handles.motif_name='Threespecies';
contents = cellstr(get(handles.growthmenu,'String'));
handles.growthmodel=contents{get(handles.growthmenu,'Value')};

handles=model_sel(handles);
guidata(hObject,handles)


% --------------------------------------------------------------------
function plotoptions_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function solutionsplot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function trajplot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function twodimplot_Callback(hObject, eventdata, handles)
set(handles.twodimplot,'Checked','on'); set(handles.threedimplot,'Checked','off')
handles.plotdim='two';
set(handles.s1_check,'value',1)
set(handles.x1_check,'value',1)
set(handles.s2_check,'value',0)
set(handles.x2_check,'value',0)
set(handles.s3_check,'value',0)
set(handles.x3_check,'value',0)
set(handles.s4_check,'value',0)
set(handles.s5_check,'value',0)
set(handles.s6_check,'value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function threedimplot_Callback(hObject, eventdata, handles)
set(handles.twodimplot,'Checked','off'); set(handles.threedimplot,'Checked','on')
handles.plotdim='three';
set(handles.s1_check,'value',1)
set(handles.x1_check,'value',1)
set(handles.s2_check,'value',1)
set(handles.x2_check,'value',0)
set(handles.s3_check,'value',0)
set(handles.x3_check,'value',0)
set(handles.s4_check,'value',0)
set(handles.s5_check,'value',0)
set(handles.s6_check,'value',0)

guidata(hObject,handles)

% --------------------------------------------------------------------
function timeplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','on'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','off')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='time';
set(handles.sol_disp,'String','Display:')
set(handles.s1_check,'Value',1,'Enable','on'); set(handles.x1_check,'Value',1,'Enable','on')
set(handles.s2_check,'Value',1,'Enable','on'); set(handles.x2_check,'Value',1,'Enable','on')
set(handles.normeig,'Visible','off','Value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function phaseplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','on'); set(handles.eigplot,'Checked','off')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='phase';
set(handles.sol_disp,'String','Select 2-3 variables')
set(handles.s1_check,'Value',0,'Enable','on'); set(handles.x1_check,'Value',0,'Enable','on')
set(handles.s2_check,'Value',0,'Enable','on'); set(handles.x2_check,'Value',0,'Enable','on')
set(handles.normeig,'Visible','off','Value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function subplot_sol_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','off'); set(handles.subplot_sol,'Checked','on')
set(handles.eigplot,'Checked','off')
handles.plotsolution='subp';
set(handles.sol_disp,'String','Display')
set(handles.s1_check,'Value',1,'Enable','on'); set(handles.x1_check,'Value',1,'Enable','on')
set(handles.s2_check,'Value',1,'Enable','on'); set(handles.x2_check,'Value',1,'Enable','on')
set(handles.normeig,'Visible','on','Value',0,'String','Overlay plots')
guidata(hObject,handles)

% --------------------------------------------------------------------
function eigplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','on')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='eig';
set(handles.sol_disp,'String','Eigenvalues')
set(handles.s1_check,'Value',0,'Enable','off'); set(handles.x1_check,'Value',0,'Enable','off')
set(handles.s2_check,'Value',0,'Enable','off'); set(handles.x2_check,'Value',0,'Enable','off')
set(handles.s3_check,'Value',0,'Enable','off'); set(handles.x3_check,'Value',0,'Enable','off')
set(handles.s4_check,'Value',0,'Enable','off'); set(handles.s5_check,'Value',0,'Enable','off');
set(handles.s6_check,'Value',0,'Enable','off');
set(handles.normeig,'Visible','on','Value',0,'String','Normalise')
guidata(hObject,handles)

% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
which_plot='sol';
plot_results;

% --- Executes on button press in lsanaly.
function lsanaly_Callback(hObject, eventdata, handles)
gvlsa=get(handles.lsanaly,'Value');
if gvlsa==0
    set(handles.lsanaly,'ForegroundColor','r')
elseif gvlsa==1
    set(handles.lsanaly,'ForegroundColor',[0.078, 0.169, 0.549])
end


% --- Executes on button press in routhcrit.
function routhcrit_Callback(hObject, eventdata, handles)
gvrhc=get(handles.routhcrit,'Value');
if gvrhc==0
    set(handles.routhcrit,'ForegroundColor','r')
elseif gvrhc==1
    set(handles.routhcrit,'ForegroundColor',[0.078, 0.169, 0.549])
end


% --- Executes on button press in jacobian_but.
function jacobian_but_Callback(hObject, eventdata, handles)
gvjac=get(handles.jacobian_but,'Value');
if gvjac==0
    set(handles.jacobian_but,'ForegroundColor','r')
elseif gvjac==1
    set(handles.jacobian_but,'ForegroundColor',[0.078, 0.169, 0.549])
end

function pert_p_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pert_p_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function sim_opts_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function sim_sp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','on'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor',[0.502 0.502 0.502])
set(handles.param1txt,'Enable','off'); set(handles.param2txt,'Enable','off')
set(handles.simparam1,'Enable','off','String','.'); set(handles.simparam2,'Enable','off','String','.'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.min_p1,'Enable','off','String','','Visible','on'); set(handles.max_p1,'Enable','off','String','','Visible','on'); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.min_p2,'Enable','off','String','','Visible','on'); set(handles.max_p2,'Enable','off','String','','Visible','on');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Enable','off','String','','Visible','on'); set(handles.step_p2,'Enable','off','String','','Visible','on','Style','edit');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Plot of trajector from initial conditions'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off')
set(handles.furth_plots,'Visible','on')
handles.simtype='single_p';

guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_mp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','on')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor','r'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Parameter 1'); set(handles.param2txt,'Enable','on','String','Parameter 2'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.simparam1,'Enable','on','String',handles.params_sim); set(handles.simparam2,'Enable','on','String',handles.params_sim); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.min_p1,'Enable','on','String','','Visible','on'); set(handles.max_p1,'Enable','on','String','','Visible','on');
set(handles.min_p2,'Enable','on','String','','Visible','on'); set(handles.max_p2,'Enable','on','String','','Visible','on');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Enable','on','String','','Visible','on'); set(handles.step_p2,'Enable','on','String','','Visible','on','Style','edit');
set(handles.min_txt,'Enable','on','Visible','on','String','Minimum'); set(handles.max_txt,'Enable','on','Visible','on','String','Maximum');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Bifurcation phase plot'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
handles.simtype='multiple_p';
set(handles.furth_plots,'Visible','off')
guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_boa_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','on'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor','b'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Variable 1'); set(handles.param2txt,'Enable','on','String','Variable 2'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.simparam1,'Enable','on','String',handles.var_names); set(handles.simparam2,'Enable','on','String',handles.var_names); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.simparam3,'Visible','off','Enable','off','String',handles.var_names);
set(handles.min_p1,'Visible','on','Enable','on','String','0'); set(handles.max_p1,'Visible','on','Enable','on','String','1');
set(handles.min_p2,'Visible','off','String',''); set(handles.max_p2,'Visible','off','String','');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Visible','on','Enable','on','String','50'); set(handles.step_p2,'Visible','off','String','','Style','edit');
set(handles.min_txt,'Visible','on','Enable','on','String','Lower Limit'); set(handles.max_txt,'Visible','on','Enable','on','String','Upper Limit');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Basin of Attraction'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.furth_plots,'Visible','off')
handles.simtype='boa';

guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_pp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','on')
set(handles.simpanel,'ForegroundColor','b'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Variable 1'); set(handles.param2txt,'Enable','on','String','Variable 2'); set(handles.var3_txt,'Visible','off','Enable','on')
set(handles.simparam1,'Enable','on','String',handles.var_names); set(handles.simparam2,'Enable','on','String',handles.var_names); set(handles.simparam3,'Visible','off','Enable','on','String','');
set(handles.min_p1,'Visible','on','Enable','on','String','0'); set(handles.max_p1,'Visible','on','Enable','on','String','1');
set(handles.min_p2,'Visible','on','Enable','on','String','0'); set(handles.max_p2,'Visible','on','Enable','on','String','1');
set(handles.min_p3,'Visible','off','Enable','on','String',''); set(handles.max_p3,'Visible','off','Enable','on','String','');
set(handles.step_p1,'Visible','on','Enable','on','String','50'); set(handles.step_p2,'Visible','on','Enable','on','String',strvcat('Fixed','Random'),'Style','popupmenu','Value',1);
set(handles.min_txt,'Visible','on','Enable','on','String','Lower IC'); set(handles.max_txt,'Visible','on','Enable','on','String','Upper IC');
set(handles.use_3v,'Visible','on','Value',0);
set(handles.uipanel6,'Title','Phase Portrait'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off')
set(handles.furth_plots,'Visible','off')
handles.simtype='pport';

guidata(hObject,handles)


% --------------------------------------------------------------------
function gui_opts_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function gui_reset_Callback(hObject, eventdata, handles)
%reset all of the values
set(handles.km1_in,'String',13);
set(handles.y1_in,'String',0.04);
set(handles.kdec1_in,'String',0.02);
set(handles.km2_in,'String',35);
set(handles.y2_in,'String',0.06);
set(handles.kdec2_in,'String',0.02);
set(handles.s1in_in,'String',5);
set(handles.ki2_in,'String',0.0000035);
set(handles.ks1_in,'String',0.3);
set(handles.ks2_in,'String',0.000025);
set(handles.n1_in,'String',3);
set(handles.n2_in,'String',2);
set(handles.d_in,'String',0.1);
set(handles.time_in,'String',1000);
set(handles.s1_init,'String',0.1);
set(handles.x1_init,'String',0.1);
set(handles.s2_init,'String',0.1);
set(handles.x2_init,'String',0.1);
set(handles.km3_in,'String',21);
set(handles.y3_in,'String',0.04);
set(handles.kdec3_in,'String',0.02);
set(handles.ks3_in,'String',0.001);
set(handles.n3_in,'String',2);
set(handles.gamma0,'String',0.43);
set(handles.gamma1,'String',0.1429);
set(handles.gamma2,'String',0.0769);
set(handles.gamma3,'String','','Visible','off')
set(handles.gamma4,'String','','Visible','off')
set(handles.gamma5,'String','','Visible','off')
set(handles.gamma6,'String','','Visible','off')

set(handles.ks32_in,'String',1e-6);
set(handles.abstol,'String',1e-8);
set(handles.reltol,'String',1e-8);
set(handles.error_val,'String','1e-6');
set(handles.step_size,'String','0.01');
set(handles.prec_num,'String','32');
set(handles.solver,'value',3);
set(handles.pert_p,'value',0.00001);
set(handles.normeig,'Visible','off','Value',0)
handles.gui_front=0;

guidata(hObject,handles)


% --------------------------------------------------------------------
function gui_restart_Callback(hObject, eventdata, handles)
close(gcbf)
MI_SIM

% --------------------------------------------------------------------
function gui_close_Callback(hObject, eventdata, handles)
close(gcbf)


% --- Executes on selection change in simparam1.
function simparam1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function simparam1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in simparam2.
function simparam2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function simparam2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p1_Callback(hObject, eventdata, handles)
minp1a = get(handles.min_p1,'String');
minp1 = str2double(minp1a);

% --- Executes during object creation, after setting all properties.
function min_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_p1_Callback(hObject, eventdata, handles)
maxp1a = get(handles.max_p1,'String');
maxp1 = str2double(maxp1a);

% --- Executes during object creation, after setting all properties.
function max_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p2_Callback(hObject, eventdata, handles)
minp2a = get(handles.min_p2,'String');
minp2 = str2double(minp2a);

% --- Executes during object creation, after setting all properties.
function min_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_p2_Callback(hObject, eventdata, handles)
maxp2a = get(handles.max_p2,'String');
maxp2 = str2double(maxp2a);
% --- Executes during object creation, after setting all properties.
function max_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function step_p1_Callback(hObject, eventdata, handles)
stepp1a = get(handles.step_p1,'String');
stepp1 = str2double(stepp1a);

% --- Executes during object creation, after setting all properties.
function step_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function step_p2_Callback(hObject, eventdata, handles)
stepp2a = get(handles.step_p2,'String');
stepp2 = str2double(stepp2a);

% --- Executes during object creation, after setting all properties.
function step_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in simparam3.
function simparam3_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function simparam3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p3_Callback(hObject, eventdata, handles)
minp3a = get(handles.min_p3,'String');
minp3 = str2double(minp3a);


% --- Executes during object creation, after setting all properties.
function min_p3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function max_p3_Callback(hObject, eventdata, handles)
maxp3a = get(handles.max_p3,'String');
maxp3 = str2double(maxp3a);

% --- Executes during object creation, after setting all properties.
function max_p3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_3v.
function use_3v_Callback(hObject, eventdata, handles)
on3d=get(handles.use_3v,'Value');

if on3d==0
    set(handles.var3_txt,'Visible','off','Enable','on')
    set(handles.simparam3,'Visible','off','Enable','on','String','');
    set(handles.min_p3,'Visible','off','Enable','on','String',''); set(handles.max_p3,'Visible','off','Enable','on','String','');
elseif on3d==1
    set(handles.var3_txt,'Visible','on','Enable','on')
    set(handles.simparam3,'Visible','on','Enable','on','String',handles.var_names);
    set(handles.min_p3,'Visible','on','Enable','on','String','0'); set(handles.max_p3,'Visible','on','Enable','on','String','1');
end
guidata(hObject,handles)

% --- Executes on button press in parallel.
function parallel_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function undock_Callback(hObject, eventdata, handles)
fig=figure;
ax=axes;
new_h=copyobj(handles.plothandle,ax,'legacy');


% --------------------------------------------------------------------
function Dimensions_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function twodphase_Callback(hObject, eventdata, handles)
if ~isempty(handles.yout_phase)
    reset(gca)
    axes(handles.trajectoryplot)
    str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
    val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
    param1=str1(val1,:); param2=str2(val2,:);
    
    h=plot(handles.yout_phase,handles.yout_phase2,handles.colr_phase,'linestyle',handles.mrk_phase); hold on
    %Plot initial conditions
    plot(handles.xlns,handles.ylns,'ks','markersize',7,'markerfacecolor','r')
    %Plot steady-states (if available), otherwise plot final values
    plot(handles.yout_phase_end,handles.yout_phase_end2,'ko','markersize',10,'markerfacecolor',handles.colre_phase)
    axis tight
    xlabel(param1); ylabel(param2)
    set(h,'linewidth',2)
    
    %Save Figure to temp_fig
    fig_name=['temp_fig/',datestr(datetime),'_phase2D_plot.pdf'];
    newfig_a=figure;
    axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
    hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
    close(newfig_a)
end

% --------------------------------------------------------------------
function threedphase_Callback(hObject, eventdata, handles)
if ~isempty(handles.yout_phase)
    reset(gca)
    axes(handles.trajectoryplot)
    str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
    val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
    param1=str1(val1,:); param2=str2(val2,:);
    
    h=plot3(handles.yout_phase,handles.yout_phase2,handles.yout_phase3,handles.colr_phase,'linestyle',handles.mrk_phase); hold on
    %Plot initial conditions
    plot3(handles.xlns,handles.ylns,handles.zlns,'ks','markersize',7,'markerfacecolor','r')
    %Plot steady-states (if available), otherwise plot final values
    plot3(handles.yout_phase_end,handles.yout_phase_end2,handles.yout_phase_end3,'ko','markersize',10,'markerfacecolor',handles.colre_phase)
    axis tight
    xlabel(param1); ylabel(param2); zlabel(str1(handles.val3,:));
    set(h,'linewidth',2); grid on
    %Save Figure to temp_fig
    fig_name=['temp_fig/',datestr(datetime),'_phase3D_plot.pdf'];
    newfig_a=figure;
    axesobject_a=copyobj(handles.trajectoryplot,newfig_a);
    hgexport(newfig_a, fig_name, hgexport('factorystyle'), 'Format', 'pdf');
    close(newfig_a)
end

% --- Executes on button press in normeig.
function normeig_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uiundock_ClickedCallback(hObject, eventdata, handles)
fig=figure;
%Check if subplot or normal figure to be copied
val_p=get(handles.subplot_sol,'Checked');
val_pp=get(handles.sim_pp,'Checked');
if strcmp(val_p,'on') && strcmp(val_pp,'off')
    kl=length(handles.plotaxes_sub);
    if kl == 3 || kl == 4
        sa = 2; sb=2;
    elseif kl == 5 || kl == 6
        sa = 3; sb=2;
    end
    for k=1:kl
        spl(k)=subplot(sa,sb,k);
        copyobj(handles.plotaxes_sub{k},spl(k),'legacy');
        ylim=get(gca,'Ylim');
        endv=handles.plotaxes_sub{1}.XData;
        set(gca,'Xlim',[0 endv(end)],'Ylim',[0 ylim(2)]);
        grid on
        legend(handles.hlegend(k).String)
        
    end
    suplabel('Time (days)');
    suplabel('Concentration (kgCOD m^{-3})','y');
else
    ax=axes;
    new_h=copyobj(handles.plothandle,ax,'legacy');
    xlabel(get(handles.xlabelhandle,'String'));
    ylabel(get(handles.ylabelhandle,'String'));
    zlabel(get(handles.zlabelhandle,'String'));
    grid on
    axis tight
end


% --------------------------------------------------------------------
function report_menu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function gen_rep_Callback(hObject, eventdata, handles)
try
    [fname,pname]=uigetfile('temp_fig/*.pdf','Select figures to include in report','MultiSelect', 'on');
    for ii=1:length(fname)
        rname{ii}=[pname,fname{ii}];
    end
    reportname=['Reports/',datestr(datetime),'_report.pdf'];
    append_pdfs(reportname,rname{:});
    set(handles.sendreport,'Enable','on')
end

% --------------------------------------------------------------------
function arch_rep_Callback(hObject, eventdata, handles)
try
    movefile('Reports/*.pdf','Reports/Archive/')
end

% --------------------------------------------------------------------
function del_rep_Callback(hObject, eventdata, handles)
%Delete all reports in Report folder
try
    delete('Reports/*.pdf')
end
%Delete all pdf images in temp_fig
try
    delete('temp_fig/*.pdf')
end

function email_add_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function email_add_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function server_add_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function server_add_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colormp.
function colormp_Callback(hObject, eventdata, handles)
axes(handles.trajectoryplot)
strcmap=get(handles.colormp,'String');
numstrcmap=get(handles.colormp,'Value');
if numstrcmap==1
    numstrcmap==2;
end
cmap = strtrim(strcmap(numstrcmap,:));
colormap(char(cmap));


% --- Executes during object creation, after setting all properties.
function colormp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function figure1_SizeChangedFcn(hObject, eventdata, handles)


% --- Executes on button press in sendreport.
function sendreport_Callback(hObject, eventdata, handles)
%Check for email and server addresses
strem=get(handles.email_add,'String');
strserv=get(handles.server_add,'String');
[fname,pname]=uigetfile('Reports/*.pdf','Select Report to send','MultiSelect', 'off');

reportfile=[pname,fname];

if ~isempty(strem) && strfind(strem,'@') && ~isempty(strserv)
    try
        setpref('Internet','E_mail',strem);
        setpref('Internet','SMTP_Server',strserv);
        sendmail(strem,'Microbial Ecology System Analysis Report',['Attached is a PDF report of the mathematical analysis of your microbial models created on ',char(datetime)],reportfile);
    end
elseif  ~isempty(strem) && strfind(strem,'@')==0 && ~isempty(strserv)
    msgbox('Enter valid e-mail address','Error: No valid e-mail');
elseif ~isempty(strem) && isempty(strserv) && strfind(strem,'@')==1 
    msgbox('Enter valid server address','Error: No valid e-mail');
else
    msgbox('Enter valid e-mail and server address','Error: No valid e-mail');
end


% --- Executes on selection change in furth_plots.
function furth_plots_Callback(hObject, eventdata, handles)
which_plot='sol';
wp=get(handles.furth_plots,'Value');
wps=get(handles.furth_plots,'String');
wp_s=strtrim(wps(wp,:));

switch wp_s
    case {'Plots...','Solutions'}
    handles.plotsolution='time';
    case {'deltaG'}
    handles.plotsolution='deltaG';
    case {'Inhibition function'}
    handles.plotsolution='inhibition';
end
plot_results;


% --- Executes during object creation, after setting all properties.
function furth_plots_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function unit_sel_Callback(hObject, eventdata, handles)
% hObject    handle to unit_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mole_sel_Callback(hObject, eventdata, handles)
% hObject    handle to mole_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cod_sel_Callback(hObject, eventdata, handles)
% hObject    handle to cod_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
