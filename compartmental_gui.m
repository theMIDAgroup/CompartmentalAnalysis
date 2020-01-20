%% @authors: Fabrice Delbary & Sara Garbarino
% contact; garbarino@dima.unige.it
function  compartmental_gui(action,varargin)
clc
warning('off')

if nargin<1,
    
    global analysis_gui
    
    if ispc==1, analysis_gui.slash='\'; else analysis_gui.slash='/'; end
    
    analysis_gui.start_path = pwd;
    analysis_gui.path = [];
    
    pathfile=[analysis_gui.start_path analysis_gui.slash 'CONFIG' analysis_gui.slash 'path.txt'];
    if exist(pathfile,'file')
        fo=fopen(pathfile,'r');
        gui_pyr.path=fgetl(fo);
        fclose(fo);
    else
        SetPath_Callback;
    end
    
    if exist(pathfile,'file')
        comp_g;
    end
    
else
    feval(action,varargin{:});
end

end


function analysis_gui = comp_g()

global analysis_gui

analysis_gui.OUTPUTfolder = [];
analysis_gui.DATAfolder = [];

analysis_gui.fig = figure('Visible','off','Position',[350,350,350,350],'MenuBar','none'); %,'Toolbar','figure');

analysis_gui.setoutput = uicontrol('Style','pushbutton','String','set OUTPUT',...
    'Position',[110,285,125,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setOUTPUT_Callback});

analysis_gui.setmodel = uicontrol('Style','popup','String','Tumor|Liver|Kidneys',...
    'Position',[110,240,86,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setMODEL_Callback});

analysis_gui.setplot = uicontrol('Style','pushbutton','String','Model',...
    'Position',[195,241.5,40,25], ...
    'HorizontalAlignment','Right',...
    'enable','off',...
    'Callback',{@setPLOT_Callback});

analysis_gui.setdata = uicontrol('Style','pushbutton','String','Choose Data',...
    'Position',[110,195,125,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setDATA_Callback});

analysis_gui.setif = uicontrol('Style','pushbutton','String','Choose IF',...
    'Position',[110,150,125,25],...
    'HorizontalAlignment','Right',...
    'Callback',{@setIF_Callback});

analysis_gui.setgut = uicontrol('Style','pushbutton','String','Choose GUT Data',...
    'Position',[110,105,125,25], ...
    'HorizontalAlignment','Right',...
    'enable','off',...
    'Callback',{@setGUT_Callback});

analysis_gui.setstart = uicontrol('Style','pushbutton','String','START',...
    'Position',[110,30,50,50],...
    'HorizontalAlignment','Right',...
    'Callback',{@setSTART_Callback});

analysis_gui.setexit = uicontrol('Style','pushbutton','String','EXIT',...
    'Position',[185,30,50,50],...
    'HorizontalAlignment','Right',...
    'Callback',{@setEXIT_Callback});

%-------------------------------------------------------------------------%
set([analysis_gui.fig,analysis_gui.setoutput,...
    analysis_gui.setmodel, analysis_gui.setplot,... 
    analysis_gui.setdata, analysis_gui.setif,analysis_gui.setgut,...
    analysis_gui.setstart, analysis_gui.setexit],...
    'Units','normalized');

set(analysis_gui.fig,'Name','Compartmetal Analysis GUI');
set(analysis_gui.fig,'NumberTitle','off');
movegui(analysis_gui.fig,'center');
%-------------------------------------------------------------------------%
set(analysis_gui.fig,'Visible','on');

end

%% Callback

function SetPath_Callback(hObject, evendata, handles)

global analysis_gui;

if isempty(analysis_gui.path),
    str='No path selected. Click Yes to procede...';
else
    str=sprintf('Current path:\n %s\n\n Click Yes to change it...',analysis_gui.path);
end

choice = questdlg(str,'path setting:', ...
    'Yes','No','Yes');
switch choice
    case 'Yes'
        path=uigetdir(analysis_gui.start_path, 'Choose path');
        
        if path == 0
        else
            analysis_gui.path=path;
            config_dir=[analysis_gui.start_path analysis_gui.slash 'CONFIG'];
            if ~exist(config_dir,'dir'), mkdir(config_dir); end
            
            fo=fopen([config_dir,analysis_gui.slash 'path.txt'],'w');
            fprintf(fo,'%s',path);
            fclose(fo);
        end
    case 'No'
end

end


function setOUTPUT_Callback(hObject, evendata, handles)

global analysis_gui

OUTPUTfolder=uigetdir(analysis_gui.start_path,'Choose OUTPUT folder');

set([analysis_gui.setoutput],'String',OUTPUTfolder,'Units','normalized');

if OUTPUTfolder == 0
else
    analysis_gui.OUTPUTfolder=OUTPUTfolder;    
    if ~exist(analysis_gui.OUTPUTfolder,'dir') 
        mkdir(analysis_gui.OUTPUTfolder);
    end
end

end


function setMODEL_Callback(hObject, evendata, handles)

global analysis_gui

aux=get(analysis_gui.setmodel,'value');

if aux == 2
    set(analysis_gui.setgut,'enable','on');    
elseif aux == 1 || aux==3
    set(analysis_gui.setgut,'enable','off');
end

set(analysis_gui.setplot,'enable','on');

end


function setPLOT_Callback(hObject, evendata, handles)

global analysis_gui

aux=get(analysis_gui.setmodel,'value');
pathfigure=[analysis_gui.start_path analysis_gui.slash 'model' analysis_gui.slash];

if aux == 1
    
    I = imread([pathfigure 'tumor.jpg']);
    figure
    imshow(I);
    
elseif aux == 2
    
    I = imread([pathfigure 'liver.jpg']);
    figure
    imshow(I);
    
elseif aux == 3
    
    I = imread([pathfigure 'kidneys.jpg']);
    figure
    imshow(I);
    
end

end


function setDATA_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_roi,location_roi]=uigetfile('*.VOISTAT','Select Data');

set([analysis_gui.setdata],'String',namefile_roi,'Units','normalized');
if namefile_roi ~= 0
    if ~iscell(namefile_roi)
        if namefile_roi==0
        else
            if isfield(analysis_gui,'namefile_roi'), analysis_gui=rmfield(analysis_gui,'namefile_roi'); end
            if isfield(analysis_gui,'location_roi'), analysis_gui=rmfield(analysis_gui,'location_roi');  end
        end
    end
    
    analysis_gui.DATA = strcat(location_roi,namefile_roi);
    aux=get(analysis_gui.setmodel,'value');
    aux_str = get(analysis_gui.setmodel,'string');
    analysis_gui.DATA = struct('model',deblank(aux_str(aux,:)));
    
    delimiter = '\t';
    startRow = 4;
    
    formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
    fileID = fopen(strcat(location_roi,namefile_roi),'r');
    dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
    fclose(fileID);
    
    analysis_gui.DATA.ROI.Study_Name = dataArray{:, 1};
    analysis_gui.DATA.ROI.Study_Name = analysis_gui.DATA.ROI.Study_Name(1,:);
    
    analysis_gui.DATA.ROI.Time_ROI = dataArray{:, 4};  
    analysis_gui.DATA.ROI.Averaged_ROI = dataArray{:, 7};
    
    analysis_gui.DATA.ROI.Volume_ROI = dataArray{:, 11};
    analysis_gui.DATA.ROI.Volume_ROI = analysis_gui.DATA.ROI.Volume_ROI(1,:);
    
else  set([analysis_gui.setdata],'String','Choose Data','Units','normalized');
end

end


function setIF_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_if,location_if]=uigetfile('*.VOISTAT','Select Data');
if namefile_if ~= 0
set([analysis_gui.setif],'String',namefile_if, ...
    'Units','normalized');

if ~iscell(namefile_if)
    if namefile_if==0
    else
        if isfield(analysis_gui,'namefile_if'), analysis_gui=rmfield(analysis_gui,'namefile_if'); end
        if isfield(analysis_gui,'location_if'), analysis_gui=rmfield(analysis_gui,'location_if');  end
    end
end

analysis_gui.IF = strcat(location_if,namefile_if);

delimiter = '\t';
startRow = 4;

formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
fileID = fopen(strcat(location_if,namefile_if),'r');
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
fclose(fileID);

analysis_gui.IF = struct('Averaged_IF',dataArray{:, 7});
analysis_gui.IF.Volume_IF = dataArray{:, 11};
analysis_gui.IF.Volume_IF = analysis_gui.IF.Volume_IF(1,:);

else  set([analysis_gui.setif],'String','Choose IF','Units','normalized');
end

end


function setGUT_Callback(hObject, evendata, handles)

global analysis_gui

[namefile_gut,location_gut]=uigetfile('*.VOISTAT','Select Data');

if namefile_gut ~= 0
set([analysis_gui.setgut],'String',namefile_gut,'Units','normalized');

if ~iscell(namefile_gut)
    if namefile_gut==0
    else
        if isfield(analysis_gui,'namefile_gut'), analysis_gui=rmfield(analysis_gui,'namefile_gut'); end
        if isfield(analysis_gui,'location_gut'), analysis_gui=rmfield(analysis_gui,'location_gut');  end
    end
end

analysis_gui.GUT = strcat(location_gut,namefile_gut);

delimiter = '\t';
startRow = 4;

formatSpec = '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%[^\n\r]';
fileID = fopen(strcat(location_gut,namefile_gut),'r');
dataArray = textscan(fileID,formatSpec,'Delimiter',delimiter,'HeaderLines',startRow-1,'ReturnOnError',false);
fclose(fileID);

analysis_gui.GUT = struct('Averaged_GUT',dataArray{:, 7});
analysis_gui.GUT.Volume_GUT = dataArray{:, 11};
analysis_gui.GUT.Volume_GUT = analysis_gui.GUT.Volume_GUT(1,:);

else 
    set([analysis_gui.setgut],'String','Choose GUT Data','Units','normalized');
end

end


function setSTART_Callback(hObject, evendata, handles)

global analysis_gui

aux = get(analysis_gui.setmodel,'value');

if aux == 1
    
    if ~exist(analysis_gui.OUTPUTfolder,'dir') 
        h = msgbox('Warning: You have to select OUTPUT folder ');
    end
    
    if ~isfield(analysis_gui,'DATA')
        h = msgbox('Warning: You have to choose Data ');
    end
    
    if  ~isfield(analysis_gui,'IF')
        hh = msgbox('Warning: You have to choose IF');
    end
    
    if exist(analysis_gui.OUTPUTfolder,'dir') && isfield(analysis_gui,'DATA') && isfield(analysis_gui,'IF') 
        tumor_gui;
    end

elseif aux == 2
    
    if ~exist(analysis_gui.OUTPUTfolder,'dir') 
        h = msgbox('Warning: You have to select OUTPUT folder ');
    end
    
    if ~isfield(analysis_gui,'DATA')
        h = msgbox('Warning: You have to choose Data ');
    end
    
    if  ~isfield(analysis_gui,'IF')
        hh = msgbox('Warning: You have to choose IF');
    end
    
    if  ~isfield(analysis_gui,'GUT')
        hhh = msgbox('Warning: You have to choose GUT Data');
    end
    
    if exist(analysis_gui.OUTPUTfolder,'dir') && isfield(analysis_gui,'DATA') && isfield(analysis_gui,'IF') && isfield(analysis_gui,'GUT')
        guts_liver_gui;
    end
     
elseif aux == 3
    
    if ~exist(analysis_gui.OUTPUTfolder,'dir') 
        h = msgbox('Warning: You have to select OUTPUT folder ');
    end
        
    if ~isfield(analysis_gui,'DATA')
        h = msgbox('Warning: You have to choose Data ');
    end

    if  ~isfield(analysis_gui,'IF')
        hh = msgbox('Warning: You have to choose IF');
    end
    
    if exist(analysis_gui.OUTPUTfolder,'dir') && isfield(analysis_gui,'DATA') && isfield(analysis_gui,'IF') 
        kidneys_gui;
    end
    
else
    msgbox('ERROR')
end

end


function setEXIT_Callback(hObject, evendata, handles)
clc
clear all
close all
% exit
end