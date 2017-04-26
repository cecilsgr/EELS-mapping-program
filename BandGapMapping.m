function varargout = BandGapMapping(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BandGapMapping_OpeningFcn, ...
                   'gui_OutputFcn',  @BandGapMapping_OutputFcn, ...
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

function BandGapMapping_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BandGapMapping (see VARARGIN)

% Choose default command line output for BandGapMapping
handles.output = hObject;

set(handles.filenametxt,'String','');

% Update handles structure
guidata(hObject, handles);

function varargout = BandGapMapping_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;



%% Panel 1: Load data

% Loading data:
function loadSIBtn_Callback(hObject, ~, handles)
t = get(handles.filenametxt,'String');
[filename, pathname] = uigetfile({'*.dm3';'*.msa';'*.txt';'*.*'},'Open file',t); 
file = sprintf('%s%s',pathname,filename);
if file ~= 0
    hbar = waitbar(0.3,'Running DM3Import.m');
    handles.data = DM3Import(file);
    waitbar(0.6,hbar,'Processing data')
    dim = size(handles.data.image_data);
    if length(dim) == 3
        handles.cube0 = reshape(handles.data.image_data,[dim(2),dim(1),dim(3)]);
        handles.cube0 = handles.data.intensity.scale*permute(handles.cube0,[2,1,3]);
        if dim(1) > dim(2)
            disp('Reshaping data')
            handles.cube0 = permute(handles.cube0,[2,1,3]);
        end
        dim = size(handles.cube0);
    end
    handles.E0 = -handles.data.zaxis.scale*handles.data.zaxis.origin:handles.data.zaxis.scale:handles.data.zaxis.scale*(dim(3)-1-handles.data.zaxis.origin);
    handles = readData(handles);
    handles.single = 0;
    setFields(hObject,handles);
    set(handles.filenametxt,'String',file);
    waitbar(1,hbar,'Done');
    close(hbar);
    updateFigs(hObject,handles);
else
    set(handles.statustxt,'String','No file');
end
guidata(hObject,handles);

function loadSpecBtn_Callback(hObject, ~, handles)
t = get(handles.filenametxt,'String');
[filename, pathname] = uigetfile({'*.dm3';'*.msa';'*.txt';'*.*'},'Open file',t); 
file = sprintf('%s%s',pathname,filename); 
if file ~= 0
    ext = filename(end-2:end);
    if strcmp(ext,'dm3')
        set(handles.filenametxt,'String',file);
        hbar = waitbar(0.3,'Reading file');
        handles.data = DM3Import(file);
        spec = handles.data.spectra_data{1,1};
        N = length(spec);
        handles.cube0 = reshape(spec,[1,1,N]).*handles.data.intensity.scale;
        handles.E0 = ((1:N)-handles.data.xaxis.origin).*handles.data.xaxis.scale;
        if max(handles.cube0(1,1,floor(N/3):N)) > 1E3
            type2 = questdlg('Normalize intensity:','Select intensity scaling','Yes','No','No');
            switch type2
                case 'Yes'; handles.cube0 = handles.cube0./max(handles.cube0);
            end
        end
    elseif strcmp(ext,'msa')
        set(handles.filenametxt,'String',file);
        hbar = waitbar(0.3,'Reading file');
        [ E, cube(1,1,:) ] = readMSA(fullfile(pathname,filename));
        handles.E0 = E';
        handles.data = cube;        % Original file
        handles.cube0 = cube;       % Original cube
    elseif strcmp(ext,'txt')
        type = questdlg('Columns in file:','Select file content','Intensity','Energy and Intensity','Intensity');
        hbar = waitbar(0.3,'Reading file');
        switch type
            case 'Intensity'
                handles.data = fscanf(fopen(file,'r'),'%f');
                cube(1,1,:) = handles.data;
                e = inputdlg({'Start energy: [eV]','End energy: [eV]'},'Figure settings',[1 50],{'0','20'});
                e1 = str2double(e{1}); 
                e2 = str2double(e{2});
                handles.E0 = e1:(e2-e1)/(length(handles.data)-1):e2;
            case 'Energy and Intensity'
                fileID = fopen(file,'r');
                handles.data = fscanf(fileID,'%f %f',[2 Inf]);
                cube(1,1,:) = handles.data(2,:);
                handles.E0 = handles.data(1,:);
        end
        if max(cube(1,1,:)) < 1 || max(cube(1,1,:)) > 1E3
            type2 = questdlg('Normalize intensity:','Select intensity scaling','Yes','No','No');
            switch type2 
                case 'Yes'; cube(1,1,:) = cube(1,1,:)./max(cube(1,1,:));
            end
        end
        handles.cube0 = cube;
    end
end

if exist('ext','var')
    waitbar(0.7,hbar,'Applying settings');
    handles.single = 1;
    handles = readData(handles);
    setFields(hObject,handles);
    set(handles.filenametxt,'String',file);
    waitbar(1,hbar,'Done'); 
    close(hbar);
    updateFigs(hObject,handles);
else
    set(handles.statustxt,'String','Error in importing file!');
end
guidata(hObject,handles);

function handles = readData(handles)
% Cubes:
handles.cube1 = handles.cube0;      % After energy align
handles.cube2 = handles.cube0;      % After background subtraction
handles.cube3 = handles.cube0;      % After binning
handles.cube4 = handles.cube0;      % After filtering
handles.cube = handles.cube0;       % After clipping
handles.E1 = handles.E0;            % After energy align
handles.E2 = handles.E0;            % After background subtraction
handles.E3 = handles.E0;            % After binning
handles.E = handles.E0;             % After clipping

function setFields(hObject,handles)
% Panel 1:
set(handles.spectrumX,'Enable','on','String',num2str(round(size(handles.cube0,2)/2)));
set(handles.spectrumY,'Enable','on','String',num2str(round(size(handles.cube0,1)/2)));
set(handles.areaX,'Enable','on','String',sprintf('1:%d',num2str(size(handles.cube0,2))));
set(handles.areaY,'Enable','on','String',sprintf('1:%d',num2str(size(handles.cube0,1))));
set(handles.viewlim1,'Enable','on','String',num2str(handles.E(1)));
set(handles.viewlim2,'Enable','on','String',num2str(handles.E(end)));
% Panel 1:
set(handles.resetBtn,'Enable','on');
set(handles.alignZLP,'Enable','on');
set(handles.applyBGsub,'Enable','on','Value',0);
set(handles.checkSubBtn,'Enable','on');
set(handles.backgroundlim1,'Enable','on')
set(handles.backgroundlim2,'Enable','on')
set(handles.binningE,'String','1','Enable','on');
set(handles.seeFiltBtn,'Enable','on');
set(handles.filterOn,'Enable','on','Value',0);
set(handles.filterN,'Enable','on');
set(handles.filterL,'Enable','on');
% Panel 3:
set(handles.singledata,'Enable','on');
set(handles.viewNfit,'Enable','on');
set(handles.direct,'Enable','on');
set(handles.fitexponent,'Enable','on');
set(handles.nvalue,'Enable','on');
set(handles.FindConstBtn,'Enable','on')
set(handles.startval,'Enable','on')
set(handles.endpointmethod,'Enable','on');
set(handles.slidingmethod,'Enable','on');
set(handles.fullmethod,'Enable','on');
set(handles.applySubFitting,'Enable','on');
% Panel 4:
set(handles.bandgaplim1,'Enable','on','String',num2str(max(handles.E(1),str2double(get(handles.bandgaplim1,'String')))));
set(handles.bandgaplim2,'Enable','on');
set(handles.fitlim1,'Enable','on','String',num2str(max(max(0,handles.E(1)),str2double(get(handles.fitlim1,'String')))));
set(handles.fitlim2,'Enable','on','String',num2str(min(handles.E(end),str2double(get(handles.fitlim2,'String')))));
set(handles.findendpt,'Enable','on');
set(handles.findslide,'Enable','on');
set(handles.endptE,'Enable','on');
set(handles.rangeE,'Enable','on');
set(handles.checkBtn,'Enable','on');
% Panel 5:
set(handles.tolerance,'Enable','on');
set(handles.r2precision,'Enable','on');
set(handles.r2counts,'Enable','on');
set(handles.i0counts,'Enable','on');
set(handles.autosave,'Enable','on');
set(handles.showfirst,'Enable','on');
set(handles.clearBtn,'Enable','on');
set(handles.FitBtn,'Enable','on');
if handles.single == 0
    set(handles.areaX,'String',sprintf('1:%d',size(handles.cube0,2)),'Enable','on');
    set(handles.areaY,'String',sprintf('1:%d',size(handles.cube0,1)),'Enable','on');
    set(handles.binningX,'String','1','Enable','on');
    set(handles.binningY,'String','1','Enable','on');
    set(handles.fulldata,'Enable','on','Value',1);
    set(handles.selecteddata,'Enable','on');
    set(handles.parallelOn,'Enable','on','Value',0);
else
    set(handles.areaX,'String','1','Enable','on');
    set(handles.areaY,'String','1','Enable','on');
    set(handles.binningX,'String','','Enable','Off'); 
    set(handles.binningY,'String','','Enable','Off'); 
    set(handles.fulldata,'Enable','off');
    set(handles.selecteddata,'Enable','off');
    set(handles.singledata,'Value',1);
    set(handles.parallelOn,'Enable','off','Value',0);
end
set(handles.resultsBtn,'Enable','off');
set(handles.workspaceBtn,'Enable','off');
set(handles.saveBtn,'Enable','off');
set(handles.statustxt,'String','File successfully imported.');
guidata(hObject,handles);

function resetAllBtn_Callback(hObject, ~, handles)
close(gcbf); 
mainGUI;

% Position:
function spectrumX_Callback(hObject, ~, handles)
x = round(str2double(get(handles.spectrumX,'String')));
if x < 1; x = 1; end
if x > size(handles.cube4,2); x = size(handles.cube4,2); end
set(hObject,'String',num2str(x));
updateFigs(hObject,handles)
guidata(hObject,handles);

function spectrumX_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function spectrumY_Callback(hObject, ~, handles)
y = round(str2double(get(handles.spectrumY,'String')));
if y < 1; y = 1; end
if y > size(handles.cube4,1); y = size(handles.cube4,1); end
set(hObject,'String',num2str(y));
updateFigs(hObject,handles)
guidata(hObject,handles);

function spectrumY_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);


% Area:
function areaX_Callback(hObject, ~, handles)
[x1,x2] = getrange(get(hObject,'String'));
if x1>x2; t=x1; x1=x2; x2=t; end
if x1<1; x1=1; end
if x2>size(handles.cube4,2); x2=size(handles.cube4,2); end
set(hObject,'String',sprintf('%d:%d',x1,x2));
guidata(hObject,handles);

function areaX_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','String','');
guidata(hObject, handles);

function areaY_Callback(hObject, ~, handles)
[y1,y2] = getrange(get(hObject,'String'));
if y1>y2; t=y1; y1=y2; y2=t; end
if y1<1; y1=1; end
if y2>size(handles.cube4,1); y2=size(handles.cube4,1); end
set(hObject,'String',sprintf('%d:%d',y1,y2));
guidata(hObject,handles);

function areaY_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','String','');
guidata(hObject, handles);


% Figures:
function updateFigs(hObject,handles)
delete(allchild(handles.axes1));
axes(handles.axes1); 
imagesc(sum(handles.cube4(:,:,:),3)./size(handles.cube4,3)); 
xlabel('x [pixels]','Fontsize',8); 
ylabel('y [pixels]','Fontsize',8); 
set(gca,'Fontsize',8,'XAxisLocation','top');
h = impoint(handles.axes1,pixel2point(handles));
setColor(h,'red')
displaySpectrum(getPosition(h),handles)
addNewPositionCallback(h,@(p) displaySpectrum(p,handles));
fcn = makeConstrainToRectFcn('impoint',get(handles.axes1,'XLim'),get(handles.axes1,'YLim'));
setPositionConstraintFcn(h,fcn);
guidata(hObject,handles);

function displaySpectrum(p,handles)
[x,y] = point2pixel(p,handles);
axes(handles.axes2);
I = squeeze(handles.cube4(y,x,:));
plot(handles.E3,I);
xlim([str2double(get(handles.viewlim1,'string')) str2double(get(handles.viewlim2,'string'))]);
xlabel('\Delta E [eV]','Fontsize',8); 
ylabel('Intensity [e^{-}]','Fontsize',8); 
set(gca,'Fontsize',8);

function point = pixel2point(handles)
x = str2double(get(handles.spectrumX,'String'));
y = str2double(get(handles.spectrumY,'String'));
point = [x y];

function [x,y] = point2pixel(p,handles)
x = round(p(1));
if x < 1; x = 1; end
if x > size(handles.cube4,2); x = size(handles.cube4,2); end
y = round(p(2));
if y < 1; y = 1; end
if y > size(handles.cube4,1); y = size(handles.cube4,1); end
set(handles.spectrumX,'String',sprintf('%d',x));
set(handles.spectrumY,'String',sprintf('%d',y));


% Energy axis:
function viewlim1_Callback(hObject, ~, handles)
lims = get(handles.axes2,'xlim');
set(handles.axes2,'xlim',[str2double(get(hObject,'string')) lims(2)]);
guidata(hObject,handles);

function viewlim1_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','String','');
guidata(hObject, handles);

function viewlim2_Callback(hObject, ~, handles)
lims = get(handles.axes2,'xlim');
set(handles.axes2,'xlim',[lims(1) str2double(get(hObject,'String'))]);
guidata(hObject,handles);

function viewlim2_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','String','');
guidata(hObject, handles);



%% Panel 2: Process data

% Reset processing:
function resetBtn_Callback(hObject, ~, handles)
handles = readData(handles);
set(handles.statustxt,'String','Processing has been reset');
set(handles.viewlim1,'String',num2str(handles.E0(1)));
set(handles.viewlim2,'String',num2str(handles.E0(end)));
updateFigs(hObject,handles);
% No ZLP aligning:
handles.shiftmap = [];
% No binning:
set(handles.binningX,'String','1');
set(handles.binningY,'String','1');
set(handles.binningE,'String','1');
% Background sub removed:
set(handles.applyBGsub,'Value',0)
% Filter off:
set(handles.filterOn,'Value',0);
% Reset data:
handles = readData(handles);
guidata(hObject,handles);

function resetBtn_CreateFcn(hObject, ~, handles)
guidata(hObject, handles);


% ZLP align:
function alignZLP_Callback(hObject, ~, handles)
cube = handles.cube0;
dim = size(cube);
E = handles.E0;
[~,n] = max(cube(:,:,1:round(dim(3)/2)),[],3);
n0 = round(mean2(n));
for i=1:size(n,1)
    for j=1:size(n,2)
        ni = n(i,j) - n0;
        if ni > 0
            I = cube(i,j,:);
            I = cat(3,I(ni+1:end),zeros(1,1,ni));
            cube(i,j,:) = I;
        elseif ni < 0
            I = cube(i,j,:);
            I = cat(3,zeros(1,1,-ni),I(1:end+ni));
            cube(i,j,:) = I;
        end
    end
end
figure()
imagesc(n.*(E(2)-E(1)));
axis image;
xlabel('Position [pixels]');
ylabel('Position [pixels]');
cb = colorbar;
ylabel(cb,'Energy shift [eV]');
title('Shift map');
a = max(max(n));
b = min(min(n));
cube = cube(:,:,a:end-b);
E = E(a:end-b) - E(n0);
handles.shiftmap = n.*(E(2)-E(1));
handles.cube1 = cube;
handles.cube2 = cube;
handles.cube3 = cube;
handles.cube4 = cube;
handles.cube = cube;
handles.E1 = E;
handles.E2 = E;
handles.E3 = E;
handles.E = E;
set(handles.statustxt,'String','ZLP aligned');
set(handles.viewlim1,'String',num2str(E(1)));
set(handles.viewlim2,'String',num2str(E(end)));
updateFigs(hObject,handles);
% Backgroud sub removed:
set(handles.applyBGsub,'Value',0)
% No binning:
set(handles.binningX,'String','1');
set(handles.binningY,'String','1');
set(handles.binningE,'String','1');
% Filter removed:
set(handles.filterOn,'Value',0);
guidata(hObject,handles);


% Manual background subtraction:
function checkSubBtn_Callback(hObject, ~, handles)
ask = {'Energy start [eV]:','Energy stop [eV]:'};
def = {get(handles.viewlim1,'String'),get(handles.viewlim2,'String')};
answer = inputdlg(ask,'Figure range',[1 50],def);
if ~isempty(answer)
    lims = [str2double(answer{1}) str2double(answer{2})];
    E1 = str2double(get(handles.backgroundlim1,'String')); 
    E2 = str2double(get(handles.backgroundlim2,'String'));
    f = figure(); subplot(2,1,1); 
    showBGfit(handles,E1,E2,lims);
    l1 = 150; l2 = 100;
    uicontrol('Style','text','String','Low range:','Position', [100 l1 120 20]);
    uicontrol('Style','text','String','High range:','Position', [100 l2 120 20]);
    handles.txt1 = uicontrol('Style','text','String',num2str(E1),'Position', [400 l1 120 20]);
    handles.txt2 = uicontrol('Style','text','String',num2str(E2),'Position', [400 l2 120 20]);
    uicontrol('Style', 'slider','Min',0,'Max',5,'Value',E1,'Position', [250 l1 120 20],'Callback', {@setBGrange1,handles,lims}); 
    uicontrol('Style', 'slider','Min',0,'Max',5,'Value',E2,'Position', [250 l2 120 20],'Callback', {@setBGrange2,handles,lims});
    uicontrol('Style', 'pushbutton', 'String', 'Copy values to GUI','Position', [400 50 130 25],'Callback', {@setBGrange,handles,f});
end
guidata(hObject, handles);

function checkSubBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function setBGrange1(hObject,~,handles,lims)
E1 = get(hObject,'Value');
E2 = str2double(get(handles.txt2,'String'));
if E1 > E2 - 0.01; E1 = E2 - 0.01; end
set(handles.txt1,'String',num2str(E1));
showBGfit(handles,E1,E2,lims);

function setBGrange2(hObject,~,handles,lims)
E2 = get(hObject,'Value');
E1 = str2double(get(handles.txt1,'String'));
if E2 < E1 + 0.01; E2 = E1 + 0.01; end
set(handles.txt2,'String',num2str(E2));
showBGfit(handles,E1,E2,lims)

function showBGfit(handles,E1,E2,lims)
cla(gca);
x = str2double(get(handles.spectrumX,'String'));
y = str2double(get(handles.spectrumY,'String'));
I = permute(handles.cube1(y,x,:),[3,1,2]);
e1 = find(handles.E1>=E1,1,'first');
e2 = find(handles.E1<=E2,1,'last');
I2 = squeeze(handles.cube1(y,x,e1:e2));
E = handles.E1(e1:e2)';
ft = fittype(@(a,r,x) a.*x.^(-r));
f = fit(E,I2,ft,'StartPoint',[1 1]);
I3 = ft(f.a,f.r,handles.E1(e1:end));
hold on;
plot(handles.E1,I,'b');
plot(handles.E1(e1:end),I3,'r');
plot([E1 E1],[0 4*I3(1)],'k',[E2,E2],[0 3*I3(1)],'k');
xlim(lims);
xlabel('Energy loss [eV]'); 
ylabel('Intensity [e^{-}]');
hold off;

function setBGrange(hObject,~,handles,fig)
set(handles.backgroundlim1,'String',get(handles.txt1,'String'));
set(handles.backgroundlim2,'String',get(handles.txt2,'String'));
set(handles.applyBGsub,'Value',0)
set(handles.statustxt,'String','Current background not subtracted.');
guidata(hObject,handles);
close(fig);

function applyBGsub_Callback(hObject, ~, handles)
if get(hObject,'Value') && strcmp(get(hObject,'Enable'),'on')
    [cube,E] = subtractBG(handles);
    set(handles.statustxt,'String','Background subtracted.');
else
    cube = handles.cube1;
    E = handles.E1;
    set(handles.statustxt,'String','Background subtraction removed.');
end
handles.cube2 = cube;
handles.cube3 = cube;
handles.cube4 = cube;
handles.cube = cube;
handles.E2 = E;
handles.E3 = E;
handles.E = E;
set(handles.binningX,'String','1');
set(handles.binningY,'String','1');
set(handles.binningE,'String','1');
set(handles.filterOn,'Value',0);
set(handles.viewlim1,'String',num2str(E(1)));
set(handles.bandgaplim1,'String',num2str(max(E(1),str2double(get(handles.bandgaplim1,'String')))));
set(handles.fitlim1,'String',num2str(max(E(1),str2double(get(handles.fitlim1,'String')))));
updateFigs(hObject,handles);
guidata(hObject,handles);

function applyBGsub_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','Value',0);
guidata(hObject, handles);

function [cube,E] = subtractBG(handles)
hbar = waitbar(0,'Applying background subtraction...');
cube = handles.cube1;
E1 = str2double(get(handles.backgroundlim1,'String'));
E2 = str2double(get(handles.backgroundlim2,'String'));
e1 = find(handles.E1>=E1,1,'first');
e2 = find(handles.E1<=E2,1,'last');
cube = zeros(size(cube,1),size(cube,2),size(cube,3)-e1+1);
Et = handles.E1(e1:e2)';
E = handles.E1(e1:end);
ft = fittype(@(a,r,x) a.*x.^(-r));
for i=1:size(cube,1)
    for j=1:size(cube,2)
        str = sprintf('Subtracting: %2.0f %%',round(((i-1)*size(cube,2)+j)/(size(cube,1)*size(cube,2))*100));
        waitbar(((i-1)*size(cube,2)+j)/(size(cube,1)*size(cube,2)),hbar,str);
        It = squeeze(handles.cube1(i,j,e1:e2));
        f = fit(Et,It,ft,'StartPoint',[1 1]);
        cube(i,j,:) = handles.cube1(i,j,e1:end) - reshape(ft(f.a,f.r,E),[1 1 length(E)]);
    end
end
waitbar(1,hbar,'Background subtracted');
close(hbar);

function backgroundlim1_Callback(hObject, ~, handles)
set(handles.applyBGsub,'Value',0)
set(handles.statustxt,'String','Current background not subtracted.');
handles.cube3 = handles.cube2;
handles.cube4 = handles.cube2;
handles.cube = handles.cube2;
handles.E3 = handles.E2;
handles.E = handles.E2;
updateFigs(hObject,handles);
guidata(hObject,handles);

function backgroundlim1_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function backgroundlim2_Callback(hObject, ~, handles)
set(handles.applyBGsub,'Value',0)
set(handles.statustxt,'String','Current background not subtracted.');
handles.cube3 = handles.cube2;
handles.cube4 = handles.cube2;
handles.cube = handles.cube2;
handles.E3 = handles.E2;
handles.E = handles.E2;
updateFigs(hObject,handles);
guidata(hObject,handles);

function backgroundlim2_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

% Binning:
function binningX_Callback(hObject, ~, handles)
a = size(handles.cube1,2);
b = round(str2double(get(hObject,'String')));
if b<1; b = 1; end
if b>a; b = a; end
set(handles.binningX,'String',num2str(b));
c = round(a/b/2);
c = c+c*(c<1);
set(handles.spectrumX,'String',num2str(c));
applybin(hObject,handles); 

function binningX_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function binningY_Callback(hObject, ~, handles)
a = size(handles.cube1,1);
b = round(str2double(get(hObject,'String')));
if b<1; b = 1; end
if b>a; b = a; end
set(handles.binningY,'String',num2str(b));
c = round(a/b/2);
c = c+c*(c<1);
set(handles.spectrumY,'String',num2str(c));
%set(handles.yVals,'String',num2str(c));
applybin(hObject,handles)

function binningY_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function binningE_Callback(hObject, ~, handles)
a = size(handles.cube1,3);
b = round(str2double(get(hObject,'String')));
if b<1; b = 1; end
if b>a; b = a; end
set(handles.binningE,'String',num2str(b));
c = round(a/b/2);
c = c+c*(c<1);
d = min(str2double(get(handles.viewlim2,'String')),c);
set(handles.viewlim2,'String',num2str(d));
applybin(hObject,handles)

function binningE_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function applybin(hObject,handles)
binX = round(str2double(get(handles.binningX,'String')));
binY = round(str2double(get(handles.binningY,'String')));
binE = round(str2double(get(handles.binningE,'String')));
cube = handles.cube2;
a = mod(size(cube,1),binY);
A = cube(1:binY:(end-a),:,:);
for i=2:binY; A = A + cube(i:binY:end-a,:,:); end
b = mod(size(cube,2),binX);
B = A(:,1:binX:end-b,:);
for i=2:binX; B = B + A(:,i:binX:end-b,:); end
c = mod(size(cube,3),binE);
C = B(:,:,1:binE:end-c);
for i=2:binE; C = C + B(:,:,i:binE:end-c); end
handles.cube3 = C;
handles.cube4 = C;
handles.cube = C;
handles.E3 = handles.E2(1:binE:end-c);
handles.E = handles.E2(1:binE:end-c);
set(handles.binningX,'String',num2str(binX));
set(handles.binningY,'String',num2str(binY));
set(handles.binningE,'String',num2str(binE));
set(handles.statustxt,'String',sprintf('Binning applied: (%d,%d,%d).',binX,binY,binE));
updateFigs(hObject,handles);
set(handles.areaX,'String',sprintf('1:%d',size(handles.cube3,2)));
set(handles.areaY,'String',sprintf('1:%d',size(handles.cube3,1)));
% Filter removed:
set(handles.filterOn,'Value',0);
guidata(hObject,handles);


% Filter:
function seeFiltBtn_Callback(hObject, ~, handles)
ask = {'Energy start [eV]:','Energy stop [eV]:'};
def = {get(handles.viewlim1,'String'),get(handles.viewlim2,'String')};
answer = inputdlg(ask,'Figure settings',[1 50],def);
if ~isempty(answer)
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    N = str2double(get(handles.filterN,'String')); 
    L = str2double(get(handles.filterL,'String'));
    I = permute(handles.cube3(y,x,:),[3,1,2]);
    fig = figure(); 
    subplot(2,1,1); plot(handles.E3,I);
    xlim([str2double(answer{1}) str2double(answer{2})]);
    xlabel('\Delta E [eV]'); 
    ylabel('Intensity [e^{-}]');
    hold on;
    showfilter(handles,N,L);
    l1 = 150; l2 = 100;
    uicontrol('Style','text','String','Maximum order:','Position', [100 l1 120 20]);
    uicontrol('Style','text','String','Length of filter:','Position', [100 l2 120 20]);
    handles.txt3 = uicontrol('Style','text','String',num2str(N),'Position', [400 l1 120 20]);
    handles.txt4 = uicontrol('Style','text','String',num2str(L),'Position', [400 l2 120 20]);
    uicontrol('Style', 'slider','Min',1,'Max',10,'Value',N,'SliderStep',[1,1]/9,'Position', [250 l1 120 20],'Callback', {@setfilter1,handles}); 
    uicontrol('Style', 'slider','Min',3,'Max',151,'Value',L,'SliderStep',[2,2]/150,'Position', [250 l2 120 20],'Callback', {@setfilter2,handles});
    uicontrol('Style', 'pushbutton', 'String', 'Copy values to GUI','Position', [400 50 130 25],'Callback', {@setfilter,handles,fig});
end
guidata(hObject, handles);

function seeFiltBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function setfilter1(hObject,~,handles)
set(handles.txt3,'String',num2str(round(get(hObject,'Value'))));
L = str2double(get(handles.txt4,'String'));
showfilter(handles,round(get(hObject,'Value')),L)

function setfilter2(hObject,~,handles)
L = round(get(hObject,'Value'));
if mod(L,2)==0; L = L+1; end
set(handles.txt4,'String',num2str(L));
N = str2double(get(handles.txt3,'String'));
showfilter(handles,N,L)

function showfilter(handles,N,L)
x = str2double(get(handles.spectrumX,'String'));
y = str2double(get(handles.spectrumY,'String'));
I = sgolayfilt(handles.cube3(y,x,:),N,L);
I = permute(I,[3,1,2]);
h = findobj(gca,'Type','line');
if length(h) > 1;   delete(h(1));   end
plot(handles.E2,I,'r');

function setfilter(hObject,~,handles,fig)
set(handles.filterN,'String',get(handles.txt3,'String'));
set(handles.filterL,'String',get(handles.txt4,'String'));
set(handles.filterOn,'Value',0)
set(handles.statustxt,'String','Current filter not applied.');
guidata(hObject,handles);
close(fig);

function filterOn_Callback(hObject, ~, handles)
if get(hObject,'Value') == 1 && strcmp(get(hObject,'Enable'),'on')
    cube = applyfilter(hObject,handles);
else
    cube = handles.cube3;
    set(handles.statustxt,'String','Filter removed.');
end
handles.cube4 = cube;
handles.cube5 = cube;
handles.cube = cube;
updateFigs(hObject,handles);
guidata(hObject,handles);

function filterOn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function cubefiltered = applyfilter(hObject,handles)
hbar = waitbar(0,'Applying filter...');
cube = handles.cube3; 
cubefiltered = zeros(size(cube));
N = round(str2double(get(handles.filterN,'String')));
L = round(str2double(get(handles.filterL,'String')));
if mod(L,2) == 0;	L = L+1; end
for i=1:size(cube,1)
    for j=1:size(cube,2)
        str = sprintf('Applying filter: %2.0f %%',round(((i-1)*size(cube,2)+j)/(size(cube,1)*size(cube,2))*100));
        waitbar(((i-1)*size(cube,2)+j)/(size(cube,1)*size(cube,2)),hbar,str);
        cubefiltered(i,j,:) = sgolayfilt(cube(i,j,:),N,L);
    end
end
set(handles.statustxt,'String','Filter applied');
set(handles.filterN,'String',num2str(N));    
set(handles.filterL,'String',num2str(L));
waitbar(1,hbar,'Filter applied');
close(hbar);
guidata(hObject,handles);

function filterN_Callback(hObject, ~, handles)
cube = applyfilter(hObject,handles);
handles.cube4 = cube;
handles.cube5 = cube;
handles.cube = cube;
updateFigs(hObject,handles);
guidata(hObject,handles);

function filterN_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function filterL_Callback(hObject, ~, handles)
cube = applyfilter(hObject,handles);
handles.cube4 = cube;
handles.cube5 = cube;
handles.cube = cube;
updateFigs(hObject,handles);
guidata(hObject,handles);

function filterL_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);



%% Panel 3: Set fitting type

% Data selection: 
function uibuttongroup2_SelectionChangedFcn(hObject, eventdata, handles)
col = [1 .8 1];
switch get(eventdata.NewValue,'Tag')
    case 'fulldata'
        set(handles.spectrumX,'BackgroundColor','white');
        set(handles.spectrumY,'BackgroundColor','white');
        set(handles.areaX,'BackgroundColor','white');
        set(handles.areaY,'BackgroundColor','white');
    case 'selecteddata'
        set(handles.spectrumX,'BackgroundColor','white');
        set(handles.spectrumY,'BackgroundColor','white');
        set(handles.areaX,'BackgroundColor',col);
        set(handles.areaY,'BackgroundColor',col);
    case 'singledata'
        set(handles.spectrumX,'BackgroundColor',col);
        set(handles.spectrumY,'BackgroundColor',col);
        set(handles.areaX,'BackgroundColor','white');
        set(handles.areaY,'BackgroundColor','white');
end
guidata(hObject, handles);

function fulldata_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function selecteddata_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function singledata_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);


% Fitting function exponent:
function viewNfit_Callback(hObject, ~, handles)
ask = {'Start exponent:','Step:','End exponent:','Show all fittings: (1/0)'};
def = {'0.15','0.05','1.5','0'};
answer = inputdlg(ask,'Testing exponent',[1 50],def);
if ~isempty(answer)
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    I = permute(handles.cube4(y,x,:),[3,1,2]);
    E = handles.E3;
    rangei = str2double(get(handles.fitlim1,'String'));
    rangef = str2double(get(handles.fitlim2,'String'));
    de = E(2)-E(1);
    i1 = find(abs(E-rangei)<de/2,1);
    i2 = find(abs(E-rangef)<de/2,1);
    inp = makefitinput(handles);
    inp.fit_method_number = 3;  % Use LLS fit
    if strcmp(answer{4},'1');  inp.show_result = 1;    end
    exponents = str2double(answer{1}):str2double(answer{2}):str2double(answer{3});
    onset = zeros(length(exponents),1);
    gof = zeros(length(exponents),1);
    lower = zeros(length(exponents),1);
    upper = zeros(length(exponents),1);
    hbar = waitbar(0,'Fitting');
    for i=1:length(exponents)
        inp.fit_exponent = exponents(i);
        fitresult1 = makefit(E(i1:i2),I(i1:i2),inp);
        onset(i) = fitresult1.onset;
        gof(i) = fitresult1.onset_gof.rsquare;
        lower(i) = onset(i) - fitresult1.onset_variation(1,1);
        upper(i) = fitresult1.onset_variation(1,2) - onset(i);
        waitbar(i/length(exponents),hbar,sprintf('%2.0f %%',round(i/length(exponents)*100)));
        if strcmp(answer{4},'1');  title(sprintf('Fit exponent: n = %1.2f',exponents(i)));  end
    end
    waitbar(1,hbar,'Done'); 
    close(hbar);
    [~,N] = max(gof);
    exponent = exponents(N);
    figure();
    subplot(2,1,1);
    hold on
    errorbar(exponents,onset,lower,upper);
    title(sprintf('Best exponent: n = %2.2f',exponent));
    xlim([exponents(1)-0.02 exponents(end)+0.02]);
    set(gca,'xticklabel',[],'box','on');
    ylabel('Fitted onset [eV]')
    subplot(2,1,2);
    plot(exponents,gof);
    xlim([exponents(1)-0.02 exponents(end)+0.02]);
    xlabel('Fit exponent{\it n}')
    ylabel('Goodness-of-fit R^2')
end
guidata(hObject, handles);

function uibuttongroup4_SelectionChangedFcn(hObject, ~, handles)
guidata(hObject, handles);

function direct_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function fitexponent_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function nvalue_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function nvalue_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) <= 0
    set(hObject,'String','1');
end
guidata(hObject,handles);


% Start value of constant:
function FindConstBtn_Callback(hObject, ~, handles)
x = str2double(get(handles.spectrumX,'String'));
y = str2double(get(handles.spectrumY,'String'));
I = permute(handles.cube4(y,x,:),[3,1,2]);
E = handles.E3;
rangei = str2double(get(handles.fitlim1,'String'));
rangef = str2double(get(handles.fitlim2,'String'));
de = E(2)-E(1);
i1 = find(abs(E-rangei)<de/2,1);
i2 = find(abs(E-rangef)<de/2,1);
inp = makefitinput(handles);
inp.fit_method_number = 3;  % Use LLS fit
inp.show_result = 1;
fitresult = makefit(E(i1:i2),I(i1:i2),inp);
title(sprintf('Fit constant: c = %f',fitresult.onset_fit.c));
guidata(hObject, handles);

function startval_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject,'String','0');
end
guidata(hObject,handles);

function startval_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);


% Fitting method:
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
col = [1 .9 0];
switch get(eventdata.NewValue,'Tag')
    case 'endpointmethod'
        set(handles.fitlim1,'BackgroundColor','white');
        set(handles.fitlim2,'BackgroundColor','white');
        set(handles.rangeE,'BackgroundColor','white');
        set(handles.endptE,'BackgroundColor',col);
        set(handles.bandgaplim1,'BackgroundColor',col);
        set(handles.bandgaplim2,'BackgroundColor',col);
        set(handles.i0counts,'Enable','on');
        set(handles.r2counts,'Enable','on');
    case 'slidingmethod'
        set(handles.fitlim1,'BackgroundColor','white');
        set(handles.fitlim2,'BackgroundColor','white');
        set(handles.rangeE,'BackgroundColor',col);
        set(handles.endptE,'BackgroundColor','white');
        set(handles.bandgaplim1,'BackgroundColor',col);
        set(handles.bandgaplim2,'BackgroundColor',col);
        set(handles.i0counts,'Enable','on');
        set(handles.r2counts,'Enable','on');
    case 'fullmethod'
        set(handles.fitlim1,'BackgroundColor',col);
        set(handles.fitlim2,'BackgroundColor',col);
        set(handles.rangeE,'BackgroundColor','white');
        set(handles.endptE,'BackgroundColor','white');
        set(handles.bandgaplim1,'BackgroundColor','white');
        set(handles.bandgaplim2,'BackgroundColor','white');
        set(handles.i0counts,'Enable','off');
        set(handles.r2counts,'Enable','off');
end
guidata(hObject, handles);

function endpointmethod_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function slidingmethod_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function fullmethod_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function applySubFitting_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function applySubFitting_Callback(hObject, ~, handles)
col = [1 .9 0];
if get(hObject,'Value')
    set(handles.bgSubrange,'String','1','Enable','on','BackgroundColor',col);
    set(handles.findBGrange,'Enable','on');
    set(handles.bgSubdist,'String','0','Enable','on','BackgroundColor',col);
    set(handles.fullmethod','Enable','off');
    if get(handles.fullmethod,'Value'); set(handles.endpointmethod,'Value',1); end
else
    set(handles.bgSubrange,'Enable','off','BackgroundColor','white');
    set(handles.findBGrange,'Enable','off');
    set(handles.bgSubdist,'Enable','off','BackgroundColor','white');
    set(handles.fullmethod','Enable','on');
end
guidata(hObject, handles);



%% Panel 4: Set fit parameters

% Onset test values:
function bandgaplim1_Callback(hObject, ~, handles)
minrange = handles.E3(3)-handles.E3(1);
if str2double(get(hObject,'String')) < handles.E3(1)
    set(hObject,'String',num2str(handles.E3(1)));
elseif str2double(get(hObject,'String')) > handles.E3(end)-minrange
    set(hObject,'String',num2str(handles.E3(end)-minrange));
end
guidata(hObject,handles);

function bandgaplim1_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[1 .9 0]);
guidata(hObject, handles);

function bandgaplim2_Callback(hObject, ~, handles)
minrange = handles.E3(3)-handles.E3(1);
if str2double(get(hObject,'String')) < handles.E3(1)+minrange
    set(hObject,'String',num2str(handles.E3(1)+minrange));
elseif str2double(get(hObject,'String')) > handles.E3(end)
    set(hObject,'String',num2str(handles.E3(end)));
end
guidata(hObject,handles);

function bandgaplim2_CreateFcn(hObject, ~, handles)
set(hObject,'BackgroundColor',[1 .9 0]);
guidata(hObject, handles);


% Fit range:
function fitlim1_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < handles.E3(1)
    set(hObject,'String',num2str(handles.E(1)));
elseif str2double(get(hObject,'String')) > handles.E3(end)
    set(hObject,'String',num2str(handles.E(1)));
end
guidata(hObject,handles);

function fitlim1_CreateFcn(hObject, ~, handles)
guidata(hObject, handles);

function fitlim2_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject,'String',num2str(handles.E3(end)));
elseif str2double(get(hObject,'String')) > handles.E3(end)
    set(hObject,'String',num2str(handles.E3(end)));
end
guidata(hObject,handles);

function fitlim2_CreateFcn(hObject, ~, handles)
guidata(hObject, handles);


% End-of-fit parameter:
function findendpt_Callback(hObject, ~, handles)
set(handles.endpointmethod,'Value',1);
E = handles.E3;
ask = {'Endpoint start [eV]:','Endpoint step [eV]:','Endpoint stop [eV]:','Show all fittings: (1/0)'};
def = {num2str(str2double(get(handles.bandgaplim2,'String'))+0.1),'0.1','5','0'};
answer = inputdlg(ask,'Testing endpoints',[1 50],def);
if ~isempty(answer)
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    I = permute(handles.cube4(y,x,:),[3,1,2]);
    inp = makefitinput(handles);
    inp.fitrange_checking = 1;
    if strcmp(answer{4},'1');  inp.show_result = 1;    end
    ept = str2double(answer{1}):str2double(answer{2}):str2double(answer{3});
    onset = zeros(length(ept),1);
    gof = zeros(length(ept),1);
    lower = zeros(length(ept),1);
    upper = zeros(length(ept),1);
    hbar = waitbar(0,'Fitting');
    for i=1:length(ept)
        inp.onset_endpoint = ept(i);
        fitresult1 = makefit(E,I,inp);
        onset(i) = fitresult1.onset;
        gof(i) = fitresult1.onset_gof.rsquare;
        lower(i) = onset(i) - fitresult1.onset_variation(1,1);
        upper(i) = fitresult1.onset_variation(1,2) - onset(i);
        waitbar(i/length(ept),hbar,sprintf('%2.0f %%',round(i/length(ept)*100)));
    end
    waitbar(1,hbar,'Done'); 
    close(hbar);
    [~,n] = max(gof.*onset);
    endpt = ept(n);
    figure();
    ax1 = subplot(2,1,1);
    e = errorbar(ept,onset,lower,upper);
    xlim([ept(1)-0.02 ept(end)+0.02]);
    ylabel('Fitted onset [eV]')
    set(e(1),'marker','x','linewidth',1);
    str = '\epsilon';
    title(sprintf('Best endpoint: %s = %2.2f eV',str,endpt));
    ax2 = subplot(2,1,2);
    plot(ept,gof);
    xlim([ept(1)-0.02 ept(end)+0.02]);
    ylabel('Goodness-of-fit R^2')
    xlabel('Endpoint \epsilon [eV]')
    linkaxes([ax1,ax2],'x')
end
guidata(hObject, handles);

function endptE_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','BackgroundColor',[1 .9 0]);
guidata(hObject, handles);

function endptE_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) > handles.E3(end)
    set(hObject,'String',num2str(handles.E3(end)));
end
guidata(hObject,handles);

function findslide_Callback(hObject, ~, handles)
set(handles.slidingmethod,'Value',1);
E = handles.E3;
ask = {'Min interval [eV]:','Step size interval [eV]:','Max interval [eV]:','Show all fittings: (1/0)'};
def = {'0.1','0.1','1.0','0'};
answer = inputdlg(ask,'Testing interval size',[1 50],def);
if ~isempty(answer)
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    I = permute(handles.cube4(y,x,:),[3,1,2]);
    inp = makefitinput(handles);
    inp.fitrange_checking = 1;
    if strcmp(answer{4},'1');  inp.show_result = 1;    end
    ranges = str2double(answer{1}):str2double(answer{2}):str2double(answer{3});
    onset = zeros(length(ranges),1);
    gof = zeros(length(ranges),1);
    lower = zeros(length(ranges),1);
    upper = zeros(length(ranges),1);
    hbar = waitbar(0,'Fitting');
    for i=1:length(ranges)
        inp.onset_interval = ranges(i);
        fitresult1 = makefit(E,I,inp);
        onset(i) = fitresult1.onset;
        gof(i) = fitresult1.onset_gof.rsquare;
        lower(i) = onset(i) - fitresult1.onset_variation(1,1);
        upper(i) = fitresult1.onset_variation(1,2) - onset(i);
        waitbar(i/length(ranges),hbar,sprintf('%2.0f %%',round(i/length(ranges)*100)));
    end
    waitbar(1,hbar,'Done'); 
    close(hbar);
    [~,n] = max(gof.*onset);
    interval = ranges(n);
    figure();
    ax1 = subplot(2,1,1);
    e = errorbar(ranges,onset,lower,upper);
    xlim([ranges(1)-0.02 ranges(end)+0.02]);
    ylabel('Fitted onset [eV]')
    set(e(1),'marker','x','linewidth',1);
    str = '\Delta';
    title(sprintf('Best interval: %s = %2.2f eV',str,interval));
    ax2 = subplot(2,1,2);
    plot(ranges,gof);
    xlim([ranges(1)-0.02 ranges(end)+0.02]);
    ylabel('Goodness-of-fit R^2')
    xlabel('Interval size \Delta [eV]')
    linkaxes([ax1,ax2],'x')
end
guidata(hObject, handles);

function rangeE_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off','BackgroundColor','white');
guidata(hObject, handles);

function rangeE_Callback(hObject, ~, handles)
minrange = handles.E3(3)-handles.E3(1);
if str2double(get(hObject,'String')) < minrange
    set(hObject,'String',num2str(minrange));
end
guidata(hObject,handles);


% Dynamic background fit range:
function findBGrange_Callback(hObject, ~, handles)
E = handles.E3;
ask = {'Background fit range length start [eV]:','Step [eV]:','Background fit range length stop [eV]:','Show all fittings: (1/0)'};
def = {'0.1','0.1','1','0'};
answer = inputdlg(ask,'Testing endpoints',[1 50],def);
if ~isempty(answer)
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    I = permute(handles.cube4(y,x,:),[3,1,2]);
    inp = makefitinput(handles);
    inp.fitrange_checking = 1;
    if strcmp(answer{4},'1');  inp.show_result = 1;    end
    BGrange = str2double(answer{1}):str2double(answer{2}):str2double(answer{3});
    onset = zeros(length(BGrange),1);
    gof = zeros(length(BGrange),1);
    lower = zeros(length(BGrange),1);
    upper = zeros(length(BGrange),1);
    hbar = waitbar(0,'Fitting');
    for i=1:length(BGrange)
        inp.background_range = BGrange(i);
        fitresult1 = makefit(E,I,inp);
        onset(i) = fitresult1.onset;
        gof(i) = fitresult1.onset_gof.rsquare;
        lower(i) = onset(i) - fitresult1.onset_variation(1,1);
        upper(i) = fitresult1.onset_variation(1,2) - onset(i);
        waitbar(i/length(BGrange),hbar,sprintf('%2.0f %%',round(i/length(BGrange)*100)));
    end
    waitbar(1,hbar,'Done'); 
    close(hbar);
    [~,n] = max(gof);
    range = BGrange(n);
    figure();
    ax1 = subplot(2,1,1);
    e = errorbar(BGrange,onset,lower,upper);
    xlim([BGrange(1)-0.02 BGrange(end)+0.02]);
    ylabel('Fitted onset [eV]')
    set(e(1),'marker','x','linewidth',1);
    str = '\xi';
    title(sprintf('Best background fit range: %s = %2.2f eV',str,range));
    ax2 = subplot(2,1,2);
    plot(BGrange,gof);
    xlim([BGrange(1)-0.02 BGrange(end)+0.02]);
    ylabel('Goodness-of-fit R^2')
    xlabel('Background fit range \xi [eV]')
    linkaxes([ax1,ax2],'x')
end
guidata(hObject, handles);

function bgSubrange_Callback(hObject, ~, handles)
minrange = handles.E3(3)-handles.E3(1);
if str2double(get(hObject,'String')) < minrange
    set(hObject,'String',num2str(minrange));
end
guidata(hObject,handles);

function bgSubrange_CreateFcn(hObject, ~, handles)
guidata(hObject, handles);

function bgSubdist_Callback(hObject, ~, handles)
minrange = handles.E3(3)-handles.E3(1);
if str2double(get(hObject,'String')) < minrange
    set(hObject,'String',num2str(minrange));
end
guidata(hObject,handles);

function bgSubdist_CreateFcn(hObject, ~, handles)
guidata(hObject, handles);


% Check parameters:
function checkBtn_Callback(hObject, eventdata, handles)
checkranges(hObject,handles);
d = dialog('Position',[300 300 250 150],'Name','Parameter checking');
uicontrol('Parent',d,'Style','text','Position',[20 80 210 40],'String','Parameter check is completed.');
uicontrol('Parent',d,'Position',[85 20 70 25],'String','Close','Callback','delete(gcf)');

function rangedeffig_Callback(hObject, ~, handles)
winopen('fitrangefig.pdf')
guidata(hObject, handles);



%% Panel 5: Run fitting

% Fit sensitivity:
function tolerance_Callback(hObject, ~, handles)

function tolerance_CreateFcn(hObject, ~, handles)

function r2precision_Callback(hObject, ~, handles)

function r2precision_CreateFcn(hObject, ~, handles)

function i0counts_Callback(hObject, ~, handles)

function i0counts_CreateFcn(hObject, ~, handles)

function r2counts_Callback(hObject, ~, handles)

function r2counts_CreateFcn(hObject, ~, handles)


% Run:
function parallelOn_Callback(hObject, ~, handles)
if get(hObject,'Value')
    set(handles.showfirst,'Enable','off','Value',0);
else
    set(handles.showfirst,'Enable','on');
end

function parallelOn_CreateFcn(hObject, ~, handles)

function autosave_Callback(hObject, ~, handles)

function showfirst_Callback(hObject, ~, handles)
if get(hObject,'Value')
    set(handles.parallelOn,'Enable','off','Value',0);
else
    set(handles.parallelOn,'Enable','on');
end

function showfirst_CreateFcn(hObject, ~, handles)

function clearBtn_Callback(hObject,~, handles)
ask = ['Are you sure you want to delete variables? ',...
    'This means binning, filtering, deconvolution and direct background ',...
    'subtraction will not be changeable.'];
answer = questdlg(ask,'Clear variables','Yes','No','No');
if strcmp(answer,'Yes')
    handles.cube0 = [];
    handles.cube1 = [];
    handles.cube2 = [];
    handles.cube3 = [];
    handles.E0 = [];
    handles.E1 = [];
    handles.E2 = []; 
    set(handles.resetBtn,'Enable','off');
    set(handles.alignZLP,'Enable','off');
    set(handles.checkSubBtn,'Enable','off');
    set(handles.applyBGsub,'Enable','off');
    set(handles.backgroundlim1,'Enable','off');
    set(handles.backgroundlim2,'Enable','off');
    set(handles.binningX,'Enable','off');
    set(handles.binningY,'Enable','off');
    set(handles.binningE,'Enable','off');
    set(handles.seeFiltBtn,'Enable','off');
    set(handles.filterOn,'Enable','off');
    set(handles.filterN,'Enable','off');
    set(handles.filterL,'Enable','off');
end
guidata(hObject,handles);

function clearBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function FitBtn_Callback(hObject, ~, handles)

% Set up fitting:
checkranges(hObject,handles);
[E,cube] = reducecube(handles);
I = transpose(reshape(cube,size(cube,1)*size(cube,2),size(cube,3)));

% Make fitinput:
fitinput = makefitinput(handles);

% Perform fit:
fitresult = cell(size(cube,1)*size(cube,2),1);
fprintf('Fitting Spectrum Image\n')

tStart = tic;
if get(handles.parallelOn,'Value')
    pool=gcp();
    parfor j=1:size(I,2)
        fitresult{j} = makefit(E,I(:,j),fitinput);
    end
    delete(pool);
else
    delete(gcp('nocreate'));
    hbar = waitbar(0,'Fitting');
    for j=1:size(I,2)
        waitbar(j/size(I,2),hbar,sprintf('%2.0f %%',round(j/size(I,2)*100)));
        fitresult{j} = makefit(E,I(:,j),fitinput);
    end
    waitbar(1,hbar,'Done'); 
    close(hbar);
end

tEnd = toc(tStart);
str = sprintf('Total fitting time: %d hours %d min %d s for %d spectra\n',floor(tEnd/60^2),round((rem(tEnd,60^2)/60)),round(rem(tEnd,60)),size(I,2));
fprintf(str);

handles.results = reshape(transpose(fitresult),size(cube,1),size(cube,2));
clear fitresult;

% Make processing struct:
if isfield(handles,'shiftmap');     process.ZLPalign = handles.shiftmap; end
if get(handles.applyBGsub,'Value'); process.subtraction = [get(handles.backgroundlim1,'String') get(handles.backgroundlim2,'String')]; end
process.binning = [get(handles.binningX,'String') get(handles.binningY,'String') get(handles.binningE,'String')];
if get(handles.filterOn,'Value');   process.filter = [get(handles.filterN,'String') get(handles.filterL,'String')]; end

% Make data struct:
handles.cube = cube;
handles.E = E;
handles.process = process;
handles.fitinput = fitinput;
handles.fittime = str;

set(handles.statustxt,'String','Fitting completed.');
set(handles.resultsBtn,'Enable','on');
set(handles.workspaceBtn,'Enable','on');
set(handles.saveBtn,'Enable','on');
guidata(hObject,handles);

if get(handles.autosave,'Value')
    dataobject = makestruct(handles);
    oldfile = dataobject.file(1:end-4);
    file = sprintf('%s.mat',oldfile);
    if file ~= 0
        save(file,'dataobject')
        disp('Saved!')
    end
    close(gcbf);
end

function FitBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);


% Additional functions:
function fitexponent = getfitexponent(handles)
if get(handles.direct,'Value'); fitexponent = 0.5;
else fitexponent = str2double(get(handle.nvalue,'String'));
end

function checkranges(hObject,handles)
fitinput = makefitinput(handles);
E_start = str2double(get(handles.fitlim1,'String'));
E_end = str2double(get(handles.fitlim2,'String'));
minrange = handles.E3(3)-handles.E3(1);

% High limit:
if floor(fitinput.fit_method_number) == 1   % Endpoint method
    if fitinput.onset_endpoint - E_end > fitinput.fit_tolerance
        E_end = fitinput.onset_endpoint;
        set(handles.fitlim2,'String',num2str(E_end));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Endpoint outside fit range',...
            'Automatically correcting high fit range limit',...
            'New value: ',E_end,' eV')));
    end
    if fitinput.testvalue_stop + minrange - fitinput.onset_endpoint > fitinput.fit_tolerance
        fitinput.testvalue_stop = fitinput.onset_endpoint - minrange;
        set(handles.bandgaplim2,'String',num2str(testvalue_stop));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Too few points fitted at high gap values',...
            'Automatically reducing band gap high limit',...
            'New value: ',fitinput.testvalue_stop,' eV')));
    end
elseif floor(fitinput.fit_method_number) == 2   % Sliding method
    if fitinput.testvalue_stop + fitinput.onset_interval - E_end > fitinput.fit_tolerance
        E_end = fitinput.testvalue_stop + fitinput.onset_interval;
        set(handles.fitlim2,'String',num2str(E_end));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Too narrow fit range: Highest band gap fit outside energy range',...
            'Automatically correcting high fit range limit',...
            'New value: ',E_end,' eV')));
    end
end

% Low limit:
if mod(fitinput.fit_method_number,1) == 0
    if E_start - fitinput.testvalue_start > fitinput.fit_tolerance
        E_start = fitinput.testvalue_start;
        set(handles.fitlim1,'String',num2str(E_start));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Too narrow range: Lowest band gap fit outside energy range',...
            'Automatically correcting low fit range limit',...
            'New value: ',E_start,' eV')));
    end
else
    if fitinput.background_range < fitinput.background_distance
        fitinput.background_range = fitinput.background_distance;
        set(handles.bgSubrange,'String',num2str(fitinput.background_range));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Background fit range larger than available distance',...
            'Automatically reducing background fit range',...
            'New value: ',fitinput.background_range,' eV')));
    end
    if E_start >  fitinput.testvalue_start - fitinput.background_range
        E_start = fitinput.testvalue_start-fitinput.background_range;
        set(handles.fitlim1,'String',num2str(E_start));
        uiwait(warndlg(sprintf('%s\n%s\n%s%1.3f%s',...
            'Lowest band gap fit background range outside energy range',...
            'Automatically reducing lower fit range',...
            'New value: ',E_start,' eV')));
    end
end
guidata(hObject,handles);

function [methodnumber,txt] = getmethod(handles)
if get(handles.endpointmethod,'Value') && ~get(handles.applySubFitting,'Value')
    methodnumber = 1;
    txt = 'Endpoint method';
elseif get(handles.endpointmethod,'Value') && get(handles.applySubFitting,'Value')
    methodnumber = 1.5;
    txt = 'Endpoint method with background subtraction';
elseif get(handles.slidingmethod,'Value') && ~get(handles.applySubFitting,'Value')
    methodnumber = 2;
    txt = 'Sliding method';
elseif get(handles.slidingmethod,'Value') && get(handles.applySubFitting,'Value')
    methodnumber = 2.5;
    txt = 'Sliding method with background subtraction';
else
    methodnumber = 3;
    txt = 'Regular least squares method';
end

function [E,cube] = reducecube(handles)
if get(handles.fulldata,'Value')
    cube = handles.cube4;
elseif get(handles.selecteddata,'Value')
    [x1,x2] = getrange(get(handles.areaX,'String'));
    [y1,y2] = getrange(get(handles.areaY,'String'));
    cube = handles.cube4(y1:y2,x1:x2,:);
else
    x = str2double(get(handles.spectrumX,'String'));
    y = str2double(get(handles.spectrumY,'String'));
    cube = handles.cube4(y,x,:);
end
rangei = str2double(get(handles.fitlim1,'String'));
rangef = str2double(get(handles.fitlim2,'String'));
de = handles.E3(2)-handles.E3(1);
startInd = find(abs(handles.E3-rangei)<de/2,1);
stopInd = find(abs(handles.E3-rangef)<de/2,1);
E = transpose(handles.E3(startInd:stopInd));
cube = cube(:,:,startInd:stopInd);

function [low, high] = getrange(str)
vals = regexp(str,'\d+(\.\d+)?','match');
if sum(size(vals))>=2
    low = str2double(vals(1));
    high = str2double(vals(end));
else
    low = 1;
    high = 1;
end

function fitinput = makefitinput(handles)
[fitmethodnumber,fitmethodtext] = getmethod(handles);
fitinput.fit_method_number = fitmethodnumber;
fitinput.fit_method = fitmethodtext;
fitinput.fit_exponent = getfitexponent(handles);
fitinput.fitrange_checking = 0;
fitinput.testvalue_start = str2double(get(handles.bandgaplim1,'String'));
fitinput.testvalue_stop = str2double(get(handles.bandgaplim2,'String'));
fitinput.fit_startvalue_const = str2double(get(handles.startval,'String'));
fitinput.onset_interval = str2double(get(handles.rangeE,'String'));
fitinput.onset_endpoint = str2double(get(handles.endptE,'String'));
fitinput.background_range = str2double(get(handles.bgSubrange,'String'));
fitinput.background_distance = str2double(get(handles.bgSubdist,'String'));
fitinput.fit_tolerance = str2double(get(handles.tolerance,'String'));
fitinput.fit_zerointensitycount = str2double(get(handles.i0counts,'String'));
fitinput.fit_decreasingR2count = str2double(get(handles.r2counts,'String'));
fitinput.fit_R2precision = str2double(get(handles.r2precision,'String'));
fitinput.show_result = get(handles.showfirst,'Value');


% Fit output:
function resultsBtn_Callback(hObject, ~, handles)
dataobject = makestruct(handles);
assignin('base','dataobjectpass',dataobject);
resultsGUI;
guidata(hObject, handles);

function resultsBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function workspaceBtn_Callback(hObject, ~, handles)
dataobject = makestruct(handles);
assignin('base','dataobject',dataobject);
guidata(hObject, handles);

function workspaceBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function saveBtn_Callback(hObject, ~, handles)
dataobject = makestruct(handles);
oldfile = dataobject.file(1:end-4);
[file,path] = uiputfile(sprintf('%s.mat',oldfile),'Save file name');
file = sprintf('%s%s',path,file);
if file ~= 0
    save(file,'dataobject')
    disp('Saved!')
end
guidata(hObject, handles);

function saveBtn_CreateFcn(hObject, ~, handles)
set(hObject,'Enable','off');
guidata(hObject, handles);

function dataobject = makestruct(handles)
dataobject.data = handles.data;
dataobject.cube = handles.cube;
dataobject.E = handles.E;
dataobject.process = handles.process;
dataobject.fitinput = handles.fitinput;
dataobject.fitresult = handles.results;
dataobject.file = get(handles.filenametxt,'String');
dataobject.fittime = handles.fittime;
dataobject.version = 'V1.0';
