function varargout = MappingResults(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MappingResults_OpeningFcn, ...
                   'gui_OutputFcn',  @MappingResults_OutputFcn, ...
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

function MappingResults_OpeningFcn(hObject, ~, handles, varargin)
if evalin('base','exist(''dataobjectpass'',''var'')')
    handles.dataobject = evalin('base','dataobjectpass');
    if strcmp(handles.dataobject.version,'V1.0')
        set(handles.resultfiletxt,'String','Results not saved');
        set(handles.origfiletxt,'String',num2str(handles.dataobject.file));
        evalin('base','clear dataobjectpass');
        handles = readData(hObject,handles);
    else
        set(handles.origfiletxt,'String','Incorrect file version!');
    end
else
    set(handles.origfiletxt,'String','No file was passed!');
end
handles.output = hObject;
guidata(hObject, handles);

function varargout = MappingResults_OutputFcn(hObject, ~, handles) 
varargout{1} = handles.output;



%% Load/show data:

function saveresults_Callback(hObject, ~, handles)
% oldfile = handles.dataobject.file(1:end-4);
% [file,path] = uiputfile(sprintf('%s.mat',oldfile),'Save file name');
% file = sprintf('%s%s',path,file);
% vari = handles.dataobject;
% save(file,'vari')
% set(handles.resultfiletxt,'String',file);
% fclose(file);
% guidata(hObject,handles);
disp()

function loadBtn_Callback(hObject, ~, handles)
[filename, pathname] = uigetfile({'*.mat';'*.*'});
file = sprintf('%s%s',pathname,filename);
if file ~= 0
    hbar = waitbar(0,'Loading file');
    fil = load(file);
    names = fieldnames(fil);
    handles.dataobject = fil.(names{1});
    waitbar(0.2,hbar,'Reading variables')
    if strcmp(handles.dataobject.version,'V3.0')
        handles.file = names;
        close(hbar);
        handles = readData(hObject,handles);
        set(handles.resultfiletxt,'String',file);
        set(handles.origfiletxt,'String',num2str(handles.dataobject.file));
    else
        error('Wrong version!')
    end
end
guidata(hObject,handles);

function handles = readData(hObject,handles)
hbar = waitbar(0.2,'Reading Variables');

% Dimensions:
N = size(handles.dataobject.cube,3);
dim = size(handles.dataobject.cube);
single = (dim(1)*dim(2))==1;

% Fit result:
onset = zeros(dim(1),dim(2));
onset_variation = zeros(dim(1),dim(2),2);
onset_const = zeros(dim(1),dim(2));
onset_gof = cell(dim(1),dim(2));
onset_test = zeros(size(handles.dataobject.cube));
onset_test_R2 = zeros(size(handles.dataobject.cube));

% Cubes:
E_subtracted = zeros(size(handles.dataobject.cube));                    % Energy for onset fitting
cube_subtracted = zeros(size(handles.dataobject.cube));                 % Vector for onset fitting
E_fitted = zeros(size(handles.dataobject.cube));                        % Energy of fitted onset
cube_fitted = zeros(size(handles.dataobject.cube));                     % Vector of fitted onset

% Check if background was fitted:
backgroundsub = mod(handles.dataobject.fitinput.fit_method_number,1)==0.5;
if backgroundsub
    cube_background_fit = zeros(size(handles.dataobject.cube));
    E_background_fit = zeros(size(handles.dataobject.cube));
    background_R2 = zeros(size(handles.dataobject.cube));
    de = handles.dataobject.E(2)-handles.dataobject.E(1);
    l = round(handles.dataobject.fitinput.background_range/de);
    background_fit = cell(dim(1),dim(2));
    background_gof = cell(dim(1),dim(2));
end

waitbar(0.4,hbar,'Transforming variables')
for i = 1:dim(1)
    for j=1:dim(2)
        onset(i,j) = handles.dataobject.fitresult{i,j}.onset;
        if isempty(handles.dataobject.fitresult{i,j}.onset_variation)
            onset_variation(i,j,:) = [0 0];
        else
            onset_variation(i,j,:) = handles.dataobject.fitresult{i,j}.onset_variation;
        end
        l1 = N - length(handles.dataobject.fitresult{i,j}.onset_test_R2);
        onset_test(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.onset_test,zeros(l1,1));
        onset_test_R2(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.onset_test_R2,zeros(l1,1));
        onset_const(i,j) = handles.dataobject.fitresult{i,j}.onset_fit.c;
        onset_gof{i,j}.sse = handles.dataobject.fitresult{i,j}.onset_gof.sse;
        onset_gof{i,j}.rsquare = handles.dataobject.fitresult{i,j}.onset_gof.rsquare;
        onset_gof{i,j}.adjrsquare = handles.dataobject.fitresult{i,j}.onset_gof.adjrsquare;
        onset_gof{i,j}.rmse = handles.dataobject.fitresult{i,j}.onset_gof.rmse;
        l3 = N - length(handles.dataobject.fitresult{i,j}.I_fitted);
        E_fitted(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.E_fitted,zeros(l3,1));
        cube_fitted(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.I_fitted,zeros(l3,1));
        l4 = N - length(handles.dataobject.fitresult{i,j}.I_subtracted);
        E_subtracted(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.E_subtracted,zeros(l4,1));
        cube_subtracted(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.I_subtracted,zeros(l4,1));
        if backgroundsub
            E_background_fit(i,j,:) = cat(1,squeeze(E_subtracted(i,j,1:l)),zeros(N-l,1));
            ft = fittype(@(a,r,x) a.*x.^(-r));
            f = handles.dataobject.fitresult{i,j}.background_fit;
            background_fit{i,j}.a = f.a;
            background_fit{i,j}.b = f.r;
            cube_background_fit(i,j,:) = ft(f.a,f.r,squeeze(E_subtracted(i,j,:)));
            background_R2(i,j,:) = cat(1,handles.dataobject.fitresult{i,j}.background_R2,zeros(l1,1));
            background_gof{i,j} = handles.dataobject.fitresult{i,j}.background_gof;
        end
    end
    waitbar(0.4+0.2*i/dim(1),hbar,'Transforming variables')
end

% Saving handles:
handles.onset = onset;
handles.onset_variation = onset_variation;
handles.onset_test = onset_test;
handles.onset_test_R2 = onset_test_R2;
handles.onset_const = onset_const;
handles.onset_gof = onset_gof;
handles.dim = dim;
handles.E = handles.dataobject.E;
handles.cuberaw = handles.dataobject.cube;
handles.backgroundsub = backgroundsub;

if backgroundsub
    handles.cube_background_fit = cube_background_fit;
    handles.E_background_fit = E_background_fit;
    handles.background_R2 = background_R2;
    handles.background_fit = background_fit;
    handles.background_gof = background_gof;
end

handles.E_subtracted = E_subtracted;
handles.cube_subtracted = cube_subtracted;
handles.E_fitted = E_fitted;
handles.cube_fitted = cube_fitted;
handles.single = single;

if single
	handles.xunits = '';
    handles.xscale = 0;
    handles.yunits = '';
    handles.yscale = 0;
else
    handles.xunits = handles.dataobject.data.xaxis.units;
    handles.xscale = handles.dataobject.data.xaxis.scale*str2double(handles.dataobject.process.binning(1));
    handles.yunits = handles.dataobject.data.yaxis.units;
    handles.yscale = handles.dataobject.data.yaxis.scale*str2double(handles.dataobject.process.binning(2));
end
waitbar(0.6,hbar,'Setting up GUI');
setFields(hObject,handles);

waitbar(0.8,hbar,'Clearing space');
handles.dataobject.data = []; 
handles.dataobject.cube = [];
handles.dataobject.fitresult = [];

waitbar(1,hbar,'Done');
close(hbar);

function setFields(hObject,handles)

% Background subtraction
if handles.backgroundsub
    x = round(size(handles.dataobject.cube,2)/2);
    y = round(size(handles.dataobject.cube,1)/2);
    set(handles.backgroundX,'String',num2str(x),'Enable','on');
    set(handles.backgroundY,'String',num2str(y),'Enable','on');
    t = handles.onset_test(handles.onset_test>0);
    set(handles.backgroundE1,'String',num2str(t(1)-handles.dataobject.fitinput.background_distance),'Enable','on');
    set(handles.backgroundE2,'String',num2str(max(handles.E)),'Enable','on');
    set(handles.showBGsub,'Enable','on');
    set(handles.showR2BG,'Enable','on');
    set(handles.showBoth,'Enable','on');
    set(handles.evalBGsub,'Enable','on');
    guidata(hObject,handles);
    updateSpectrumBG(handles);
else
    set(handles.backgroundX,'String','','Enable','off');
    set(handles.backgroundY,'String','','Enable','off');
    set(handles.backgroundE1,'String','','Enable','off');
    set(handles.backgroundE2,'String','','Enable','off');
    set(handles.showBGsub,'Enable','off');
    set(handles.showR2BG,'Enable','off');
    set(handles.showBoth,'Enable','off');
    set(handles.evalBGsub,'Enable','off');
    set(handles.axes3,'box','on','xtick',[],'ytick',[])
end

% Pixel fields
set(handles.onsetX,'String',num2str(round(size(handles.dataobject.cube,2)/2)));
set(handles.onsetY,'String',num2str(round(size(handles.dataobject.cube,1)/2)));
set(handles.onsetE1,'String',num2str(handles.onset_test(1,1,1)));
set(handles.onsetE2,'String',num2str(max(handles.E)));
guidata(hObject,handles);
updateSpectrum(handles);

% Area fields
if ~handles.single
    set(handles.statisticsX1,'String','1');
    set(handles.statisticsX2,'String',num2str(size(handles.dataobject.cube,2)));
    set(handles.statisticsY1,'String','1');
    set(handles.statisticsY2,'String',num2str(size(handles.dataobject.cube,1)));
    set(handles.statisticsE1,'String',num2str(min(min(handles.onset))));
    set(handles.statisticsE2,'String',num2str(max(max(handles.onset))));
    set(handles.showMap,'Enable','on');
    set(handles.gapStats,'Enable','on');
    set(handles.xProfile,'Enable','on');
    set(handles.yProfile,'Enable','on');
    set(handles.showLine,'Enable','on');
    set(handles.histBins,'Enable','on');
    set(handles.gaussianHist,'Enable','on');
    set(handles.showHistogram,'Enable','on');
    set(handles.getHistogram,'Enable','on');
    guidata(hObject,handles);
    updateFig(hObject,handles);
else
	set(handles.statisticsX1,'String','','Enable','off');
    set(handles.statisticsX2,'String','','Enable','off');
    set(handles.statisticsY1,'String','','Enable','off');
    set(handles.statisticsY2,'String','','Enable','off');
    set(handles.statisticsE1,'String','','Enable','off');
    set(handles.statisticsE2,'String','','Enable','off');
    set(handles.showMap,'Enable','off');
    set(handles.gapStats,'Enable','off');
    set(handles.xProfile,'Enable','off');
    set(handles.yProfile,'Enable','off');
    set(handles.showLine,'Enable','off');
    set(handles.histBins,'Enable','off');
    set(handles.gaussianHist,'Enable','off');
    set(handles.showHistogram,'Enable','off');
    set(handles.getHistogram,'Enable','off');
    set(handles.axes1,'box','on','xtick',[],'ytick',[])
end
guidata(hObject,handles);



%% Figure settings:

function fontsize_Callback(hObject, ~, handles)

function fontsize_CreateFcn(hObject, ~, handles)
set(hObject,'String','10');

function cmap_Callback(hObject, ~, handles)
axes(handles.axes1); 
colormap(get(handles.cmap,'String'));

function cmap_CreateFcn(hObject, ~, handles)
set(hObject,'String','Jet');



%% Background subtraction:

% Spectrum:
function updateSpectrumBG(handles)
axes(handles.axes3);
cla();
[xp,yp] = getpixelBG(handles);
e1 = handles.E;
E1 = e1(e1~=0);
I1 = squeeze(handles.cuberaw(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_background_fit(yp,xp,:));
E2 = e2(e2~=0);
i2 = squeeze(handles.cube_background_fit(yp,xp,:));
I2 = i2(e2~=0);
e3 = squeeze(handles.E_subtracted(yp,xp,:));
E3 = e3(e3~=0);
I2r = i2(e3~=0);
I3 = squeeze(handles.cube_subtracted(yp,xp,:));
I3 = I3(e3~=0);
hold on
a = area(E2,I2);
plot(E1,I1,E3,I2r,E3,I3);
hold off
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.backgroundE2,'String'))]);
xlabel('Energy loss [eV]','Fontsize',8); 
ylabel('Intensity [a.u.]','Fontsize',8); 
set(a,'FaceColor',[1 0.5 0.4],'EdgeColor','none');

function [x,y] = getpixelBG(handles)
x = str2double(get(handles.backgroundX,'String'));
y = str2double(get(handles.backgroundY,'String'));

% Position:
function backgroundX_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 1
    set(hObject,'String','1');
elseif str2double(get(hObject,'String')) > handles.dim(2)
    set(hObject,'String',num2str(handles.dim(2)));
end
updateSpectrumBG(handles)

function backgroundX_CreateFcn(hObject, ~, handles)
set(hObject,'String','');

function backgroundY_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 1
    set(hObject,'String','1');
elseif str2double(get(hObject,'String')) > handles.dim(1)
    set(hObject,'String',num2str(handles.dim(1)));
end
updateSpectrumBG(handles)

function backgroundY_CreateFcn(hObject, ~, handles)
set(hObject,'String','');

% Energy range:
function backgroundE1_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < min(handles.E)
    set(hObject,'String',num2str(min(handles.E)));
end
axes(handles.axes3);
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.backgroundE2,'String'))]); 

function backgroundE1_CreateFcn(hObject, ~, handles)
set(hObject,'String','');

function backgroundE2_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) > max(handles.E)
    set(hObject,'String',num2str(max(handles.E)));
end
axes(handles.axes3);
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.backgroundE2,'String'))]); 

function backgroundE2_CreateFcn(hObject, ~, handles)
set(hObject,'String','');

% Figures:
function showBGsub_Callback(hObject, ~, handles)
[xp,yp] = getpixelBG(handles);
e1 = handles.E;
E1 = e1(e1~=0);
I1 = squeeze(handles.cuberaw(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_background_fit(yp,xp,:));
E2 = e2(e2~=0);
i2 = squeeze(handles.cube_background_fit(yp,xp,:));
I2 = i2(e2~=0);
e3 = squeeze(handles.E_subtracted(yp,xp,:));
E3 = e3(e3~=0);
I2r = i2(e3~=0);
I3 = squeeze(handles.cube_subtracted(yp,xp,:));
I3 = I3(e3~=0);
figure();
subplot(1,2,1);
hold on
a = area(E2,I2);
p = plot(E1,I1,E3,I2r,E3,I3);
hold off
set(a,'FaceColor',[1 0.5 0.4],'EdgeColor','none');
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.backgroundE2,'String'))]);
xlabel('Energy loss [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('Intensity [e^{-}]','Fontsize',str2double(get(handles.fontsize,'String'))); 
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 
l = legend(p,{'Spectrum','Background','Inelastic'});
set(l,'Location','northwest','Box','off')
gof = handles.background_gof{yp,xp};

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 350 200 20],'fontweight','bold',...
    'String','Fitting output:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 280 110 70],'String',sprintf(...
    'Model: \nFit constant: \nExponent: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 280 90 70],'String',sprintf(...
    ' ''Power law''\n %1.3f a.u.\n %1.3f a.u.',...
    handles.background_fit{yp,xp}.a,handles.background_fit{yp,xp}.b));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 250 200 20],'fontweight','bold',...
    'String','Goodness of fit:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 180 110 70],'String',sprintf(...
    'Sum square error: \nR^2: \nAdjusted R^2: \nRMS Error: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 180 90 70],'String',sprintf(...
    ' %2.2f\n %1.3f\n %1.3f\n %1.3f',...
    gof.sse,gof.rsquare,gof.adjrsquare,gof.rmse));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 150 200 20],'fontweight','bold',...
    'String','Position:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 80 110 70],'String',sprintf(...
    'x position: \nx distance: \ny position: \ny distance:'));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 80 90 70],'String',sprintf(...
    ' %d px\n %3.1f %s\n %d px\n %3.1f %s',...
    xp,xp*handles.xscale,handles.xunits,yp,yp*handles.yscale,handles.yunits));

function showR2BG_Callback(hObject, ~, handles)
[xp,yp] = getpixelBG(handles);
figure();
gap = squeeze(handles.onset_test(yp,xp,:));
background_R2 = squeeze(handles.background_R2(yp,xp,:));
R2 = squeeze(handles.onset_test_R2(yp,xp,:));
gaps = gap(gap~=0);
background_R2 = background_R2(gap~=0);
R2 = R2(gap~=0);
plot(gaps,background_R2,gaps,R2,handles.onset(yp,xp),handles.background_gof{yp,xp}.rsquare,'k*',handles.onset(yp,xp),handles.onset_gof{yp,xp}.rsquare,'k*');
xlabel('Test band gap [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('R^2','Fontsize',str2double(get(handles.fontsize,'String'))); 
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 
l = legend('Background fit','Band gap fit');
set(l,'Location','best','Box','off')

function showBoth_Callback(hObject, ~, handles)
figure();
[xp,yp] = getpixelBG(handles);
e1 = handles.E;
E1 = e1(e1~=0);
I1 = squeeze(handles.cuberaw(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_background_fit(yp,xp,:));
E2 = e2(e2~=0);
i2 = squeeze(handles.cube_background_fit(yp,xp,:));
I2 = i2(e2~=0);
e3 = squeeze(handles.E_subtracted(yp,xp,:));
E3 = e3(e3~=0);
I2r = i2(e3~=0);
I3 = squeeze(handles.cube_subtracted(yp,xp,:));
I3 = I3(e3~=0);
e4 = squeeze(handles.E_fitted(yp,xp,:));
E4 = e4(e4~=0);
I4 = squeeze(handles.cube_fitted(yp,xp,:));
I4 = I4(e4~=0);
hold on
a = area(E2,I2);
p = plot(E1,I1,E3,I2r,E3,I3,E4,I4);
plot(E1,zeros(length(E1),1),'k--');
hold off
set(a,'FaceColor',[1 0.5 0.4],'EdgeColor','none');
xlabel('Energy loss [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('Intensity [e^{-}]','Fontsize',str2double(get(handles.fontsize,'String'))); 
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.onsetE2,'String'))]);
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 
title(sprintf('Band gap: %2.2f eV',handles.onset(yp,xp)));
l = legend(p,'Spectrum','Background','Inelastic','Fitting');
set(l,'Location','northwest','Box','off')

function evalBGsub_Callback(hObject, ~, handles)
fig = figure(); 
fig1 = subplot(2,2,1);
if sum(size(handles.onset))>1
    imagesc(handles.onset,[min(min(handles.onset)) max(max(handles.onset))]); 
    colormap(get(handles.cmap,'String')); 
    axis image;  
    xlabel('x [pixels]'); ylabel('y [pixels]')
end
pos = get(fig,'Position'); 
set(fig,'Position',[pos(1) pos(2) 1.5*pos(3) pos(4)]);
h = impoint(fig1,[1 1]);
fcn = makeConstrainToRectFcn('impoint',get(fig1,'XLim'),get(fig1,'YLim'));
setPositionConstraintFcn(h,fcn); 
DisplayFittingWinBG(getPosition(h),handles)
addNewPositionCallback(h,@(p)DisplayFittingWinBG(p,handles));

function DisplayFittingWinBG(point,handles)
xp = round(point(1,1));  
yp = round(point(1,2));
subplot(1,2,2);
cla();
e1 = handles.E;
E1 = e1(e1~=0);
I1 = squeeze(handles.cuberaw(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_background_fit(yp,xp,:));
E2 = e2(e2~=0);
i2 = squeeze(handles.cube_background_fit(yp,xp,:));
I2 = i2(e2~=0);
e3 = squeeze(handles.E_subtracted(yp,xp,:));
E3 = e3(e3~=0);
I2r = i2(e3~=0);
I3 = squeeze(handles.cube_subtracted(yp,xp,:));
I3 = I3(e3~=0);
hold on
a = area(E2,I2);
p = plot(E1,I1,E3,I2r,E3,I3);
hold off
xlim([str2double(get(handles.backgroundE1,'String')) str2double(get(handles.backgroundE2,'String'))]);
xlabel('Energy loss [eV]')
ylabel('Intensity [a.u.]')
l = legend(p,'Spectrum','Background','Inelastic');
set(l,'location','best','box','off');
set(a,'FaceColor',[1 0.5 0.4],'EdgeColor','none');
set(l,'Location','northwest','Box','off')
gof = handles.background_gof{yp,xp};

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 190 200 20],'fontweight','bold',...
    'String','Fitting output:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 110 110 70],'String',sprintf(...
    'Model: \nFit constant: \nExponent: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[280 110 90 70],'String',sprintf(...
    ' ''Power law''\n %1.3f a.u.\n %1.3f a.u.',...
    handles.background_fit{yp,xp}.a,handles.background_fit{yp,xp}.b));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 80 200 20],'fontweight','bold',...
    'String','Goodness of fit:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 10 110 70],'String',sprintf(...
    'Sum square error: \nR^2: \nAdjusted R^2: \nRMS Error: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[280 10 90 70],'String',sprintf(...
    ' %2.2f\n %1.3f\n %1.3f\n %1.3f',...
    gof.sse,gof.rsquare,gof.adjrsquare,gof.rmse));



%% Fit evaluation:

% Spectrum:
function updateSpectrum(handles)
axes(handles.axes2);
[xp,yp] = getpixel(handles);
e1 = squeeze(handles.E_subtracted(yp,xp,:));
E1 = e1(e1~=0);
I1 = squeeze(handles.cube_subtracted(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_fitted(yp,xp,:));
E2 = e2(e2~=0);
I2 = squeeze(handles.cube_fitted(yp,xp,:));
I2 = I2(e2~=0);
plot(E1,I1,E2,I2);
xlabel('Energy loss [eV]','Fontsize',8); 
ylabel('Intensity [a.u.]','Fontsize',8);
set(gca,'Fontsize',8);
xlim([str2double(get(handles.onsetE1,'String')) str2double(get(handles.onsetE2,'String'))]);

function [x,y] = getpixel(handles)
x = str2double(get(handles.onsetX,'String'));
y = str2double(get(handles.onsetY,'String'));


% Position:
function onsetX_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 1
    set(hObject,'String','1');
elseif str2double(get(hObject,'String')) > handles.dim(2)
    set(hObject,'String',num2str(handles.dim(2)));
end
updateSpectrum(handles)

function onsetX_CreateFcn(hObject, ~, handles)

function onsetY_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < 1
    set(hObject,'String','1');
elseif str2double(get(hObject,'String')) > handles.dim(1)
    set(hObject,'String',num2str(handles.dim(1)));
end
updateSpectrum(handles)

function onsetY_CreateFcn(hObject, ~, handles)


% Energy range:
function onsetE1_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) < min(handles.E)
    set(hObject,'String',num2str(min(handles.E)));
end
axes(handles.axes2);
xlim([str2double(get(handles.onsetE1,'String')) str2double(get(handles.onsetE2,'String'))]); 

function onsetE1_CreateFcn(hObject, ~, handles)

function onsetE2_Callback(hObject, ~, handles)
if str2double(get(hObject,'String')) > max(handles.E)
    set(hObject,'String',num2str(max(handles.E)));
end
axes(handles.axes2);
xlim([str2double(get(handles.onsetE1,'String')) str2double(get(handles.onsetE2,'String'))]); 

function onsetE2_CreateFcn(hObject, ~, handles)


% Figures: 
function showspectrum_Callback(hObject, ~, handles)
[xp,yp] = getpixel(handles);
e1 = squeeze(handles.E_subtracted(yp,xp,:));
E1 = e1(e1~=0);
I1 = squeeze(handles.cube_subtracted(yp,xp,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_fitted(yp,xp,:));
E2 = e2(e2~=0);
I2 = squeeze(handles.cube_fitted(yp,xp,:));
I2 = I2(e2~=0);
figure();
p = plot(E1,I1,'b',E2,I2,'r',E1,zeros(length(E1),1),'k--');
set(p(2),'linewidth',2)
if handles.onset_variation(yp,xp,1)==0&&handles.onset_variation(yp,xp,2)==0
    title(sprintf('E_{onset} = %2.2f eV',handles.onset(yp,xp)));
else
    title(sprintf('E_{onset} = %2.2f eV (%2.2f eV, %2.2f eV)',handles.onset(yp,xp),handles.onset_variation(yp,xp,:)));
end
xlim([str2double(get(handles.onsetE1,'String')) str2double(get(handles.onsetE2,'String'))]);
xlabel('Energy loss [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('Intensity [e^{-}]','Fontsize',str2double(get(handles.fontsize,'String'))); 
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 
l = legend({'Spectrum','Fitting'});
set(l,'Location','northwest','Box','off')

function showR2_Callback(hObject, ~, handles)
[xp,yp] = getpixel(handles);
figure();
gap = squeeze(handles.onset_test(yp,xp,:));
R2 = squeeze(handles.onset_test_R2(yp,xp,:));
ci = squeeze(handles.onset_variation(yp,xp,:));
gaps = gap(gap~=0);
R2 = R2(gap~=0);
plot(gaps,R2,handles.onset(yp,xp),handles.onset_gof{yp,xp}.rsquare,'*');
if ci(1)>0 && ci(2)>0
    y = handles.onset_gof{yp,xp}.rsquare-handles.dataobject.fitinput.fit_R2precision;
    hold on;
    plot([gaps(1) gaps(end)],[y y],'k--',ci,[y,y],'k*',...
        [ci(1) ci(1)],[y 0],'k--',[ci(2) ci(2)],[y 0],'k--');
    title(sprintf('E_{onset} = %2.2f eV (%2.2f eV, %2.2f eV)',handles.onset(yp,xp),ci));
else
    title(sprintf('E_{onset} = %2.2f eV',handles.onset(yp,xp)));
end
xlabel('Test onset E_g'' [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('R^2','Fontsize',str2double(get(handles.fontsize,'String'))); 
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 

function displayStats_Callback(hObject, ~, handles)
[x,y] = getpixel(handles);
e1 = squeeze(handles.E_subtracted(y,x,:));
E1 = e1(e1~=0);
I1 = squeeze(handles.cube_subtracted(y,x,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_fitted(y,x,:));
E2 = e2(e2~=0);
I2 = squeeze(handles.cube_fitted(y,x,:));
I2 = I2(e2~=0);
figure();
subplot(1,2,1);
plot(E1,I1,E2,I2);
xlim([str2double(get(handles.onsetE1,'String')) str2double(get(handles.onsetE2,'String'))]);
xlabel('Energy loss [eV]','Fontsize',str2double(get(handles.fontsize,'String'))); 
ylabel('Intensity [e^{-}]','Fontsize',str2double(get(handles.fontsize,'String'))); 
set(gca,'Fontsize',str2double(get(handles.fontsize,'String'))); 
onset_gof = handles.onset_gof{y,x};

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 350 200 20],'fontweight','bold',...
    'String','Fitting output:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 280 110 70],'String',sprintf(...
    'Fit constant: \nBand onset: \nLow conf. bound: \nHigh conf. bound:'));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 280 90 70],'String',sprintf(...
    ' %1.3f a.u.\n %1.3f eV\n %1.3f eV\n %1.3f eV',...
    handles.onset_const(y,x),handles.onset(y,x),squeeze(handles.onset_variation(y,x,:)) ));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 250 200 20],'fontweight','bold',...
    'String','Goodness of fit:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 180 110 70],'String',sprintf(...
    'Sum square error: \nR^2: \nAdjusted R^2: \nRMS Error: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 180 90 70],'String',sprintf(...
    ' %2.2f\n %1.3f\n %1.3f\n %1.3f',...
    onset_gof.sse,onset_gof.rsquare,onset_gof.adjrsquare,onset_gof.rmse));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 150 200 20],'fontweight','bold',...
    'String','Position:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[300 80 110 70],'String',sprintf(...
    'x position: \nx distance: \ny position: \ny distance:'));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[410 80 90 70],'String',sprintf(...
    ' %d px\n %3.1f %s\n %d px\n %3.1f %s',...
    x,x*handles.xscale,handles.xunits,y,y*handles.yscale,handles.yunits));

function fitStats_Callback(hObject, ~, handles)
fig = figure(); 
fig1 = subplot(2,2,1);
if sum(size(handles.onset))>1
    imagesc(handles.onset,[min(min(handles.onset)) max(max(handles.onset))]); 
    colormap(get(handles.cmap,'String')); 
    axis image;  
    xlabel('x [pixels]'); ylabel('y [pixels]')
end
pos = get(fig,'Position'); 
set(fig,'Position',[pos(1) pos(2) 1.5*pos(3) pos(4)]);
h = impoint(fig1,[1 1]);
fcn = makeConstrainToRectFcn('impoint',get(fig1,'XLim'),get(fig1,'YLim'));
setPositionConstraintFcn(h,fcn); 
DisplayFittingWin(getPosition(h),handles)
addNewPositionCallback(h,@(p)DisplayFittingWin(p,handles));

function DisplayFittingWin(point,handles)
x = round(point(1,1));  
y = round(point(1,2));
subplot(1,2,2);
cla();
e1 = squeeze(handles.E_subtracted(y,x,:));
E1 = e1(e1~=0);
I1 = squeeze(handles.cube_subtracted(y,x,:));
I1 = I1(e1~=0);
e2 = squeeze(handles.E_fitted(y,x,:));
E2 = e2(e2~=0);
I2 = squeeze(handles.cube_fitted(y,x,:));
I2 = I2(e2~=0);
plot(E1,I1,E2,I2,E1,zeros(length(E1),1),'k--');
xlim([min(min(handles.onset))-0.2 E1(end)])
xlabel('Energy [eV]')
ylabel('Intensity [e^{-}]')
onset_gof = handles.onset_gof{y,x};

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 190 200 20],'fontweight','bold',...
    'String','Fitting output:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 110 110 70],'String',sprintf(...
    'Fit constant: \nBand onset: \nLow conf. bound: \nHigh conf. bound:'));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[280 110 90 70],'String',sprintf(...
    ' %1.3f a.u.\n %1.3f eV\n %1.3f eV\n %1.3f eV',...
    handles.onset_const(y,x),handles.onset(y,x),squeeze(handles.onset_variation(y,x,:)) ));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 80 200 20],'fontweight','bold',...
    'String','Goodness of fit:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[150 10 110 70],'String',sprintf(...
    'Sum square error: \nR^2: \nAdjusted R^2: \nRMS Error: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[280 10 90 70],'String',sprintf(...
    ' %2.2f\n %1.3f\n %1.3f\n %1.3f',...
    onset_gof.sse,onset_gof.rsquare,onset_gof.adjrsquare,onset_gof.rmse));



%% Band gap statistics:

% Map:
function updateFig(hObject,handles)
delete(allchild(handles.axes1));
axes(handles.axes1); 
lim1 = str2double(get(handles.statisticsE1,'String'));
lim2 = str2double(get(handles.statisticsE2,'String'));
if lim2 <= lim1 % If only 1 onset value
    lim1 = lim1-0.25; set(handles.statisticsE1,'String',num2str(lim1));
    lim2 = lim1+0.5; set(handles.statisticsE2,'String',num2str(lim2));
end
imagesc(handles.onset,[lim1 lim2]);
colormap(get(handles.cmap,'String'));
xlabel('x [pixels]','Fontsize',8); 
ylabel('y [pixels]','Fontsize',8); 
set(gca,'Fontsize',8,'XAxisLocation', 'top')
guidata(hObject,handles);

function [x1,x2,y1,y2] = getarea(handles)
x1 = str2double(get(handles.statisticsX1,'String'));
x2 = str2double(get(handles.statisticsX2,'String'));
y1 = str2double(get(handles.statisticsY1,'String'));
y2 = str2double(get(handles.statisticsY2,'String'));

function impos = coord2impos(handles)
[x1,x2,y1,y2] = getarea(handles);
impos = [x1-0.5 y1-0.5 x2-x1+1 y2-y1+1];


% Area:
function statisticsX1_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if val < 1
    set(hObject,'String','1');
elseif val > handles.dim(2)
    set(hObject,'String',num2str(handles.dim(2)));
end
if val > str2double(get(handles.statisticsX2,'String'))
    set(hObject,'String',get(handles.statisticsX2,'String'));
    set(handles.statisticsX2,'String',num2str(val));
end
updateFig(hObject,handles)

function statisticsX1_CreateFcn(hObject, ~, handles)

function statisticsX2_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if val < 1
    set(hObject,'String','1');
elseif val > handles.dim(2)
    set(hObject,'String',num2str(handles.dim(2)));
end
if val < str2double(get(handles.statisticsX1,'String'))
    set(hObject,'String',get(handles.statisticsX1,'String'));
    set(handles.statisticsX1,'String',num2str(val));
end
updateFig(hObject,handles)

function statisticsX2_CreateFcn(hObject, ~, handles)

function statisticsY1_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if val < 1
    set(hObject,'String','1');
elseif val > handles.dim(1)
    set(hObject,'String',num2str(handles.dim(1)));
end
if val > str2double(get(handles.statisticsY2,'String'))
    set(hObject,'String',get(handles.statisticsY2,'String'));
    set(handles.statisticsY2,'String',num2str(val));
end
updateFig(hObject,handles)

function statisticsY1_CreateFcn(hObject, ~, handles)

function statisticsY2_Callback(hObject, ~, handles)
val = str2double(get(hObject,'String'));
if val < 1
    set(hObject,'String','1');
elseif val > handles.dim(1)
    set(hObject,'String',num2str(handles.dim(1)));
end
if val < str2double(get(handles.statisticsY1,'String'))
    set(hObject,'String',get(handles.statisticsY1,'String'));
    set(handles.statisticsY1,'String',num2str(val));
end
updateFig(hObject,handles)

function statisticsY2_CreateFcn(hObject, ~, handles)


% Band gap limits:
function statisticsE1_Callback(hObject, ~, handles)
axes(handles.axes1);
caxis([str2double(get(handles.statisticsE1,'String')) str2double(get(handles.statisticsE2,'String'))]); 

function statisticsE1_CreateFcn(hObject, ~, handles)

function statisticsE2_Callback(hObject, ~, handles)
axes(handles.axes1);
caxis([str2double(get(handles.statisticsE1,'String')) str2double(get(handles.statisticsE2,'String'))]); 

function statisticsE2_CreateFcn(hObject, ~, handles)


% Line profile:
function xProfile_Callback(hObject, ~, handles)
if get(hObject,'Value'); 
    set(handles.yProfile,'Value',0);
end
guidata(hObject,handles);

function yProfile_Callback(hObject, ~, handles)
if get(hObject,'Value'); 
    set(handles.xProfile,'Value',0);
end
guidata(hObject,handles);

function errorbarOn_Callback(hObject, ~, handles)

function errorStd_Callback(hObject, ~, handles)
if get(hObject,'Value'); 
    set(handles.pointCI,'Value',0);
end
guidata(hObject,handles);

function pointCI_Callback(hObject, ~, handles)
if get(hObject,'Value'); 
    set(handles.errorStd,'Value',0);
end
guidata(hObject,handles);

function showLine_Callback(hObject, ~, handles)
[x1,x2,y1,y2] = getarea(handles);
if get(handles.xProfile,'Value')
    if get(handles.pointCI,'Value')
        y1 = round((y1+y2)/2);
        gapvec = handles.onset(y1,x1:x2);
        errvec(:,1) = squeeze(handles.onset_variation(y1,x1:x2,1))-gapvec;
        errvec(:,2) = squeeze(handles.onset_variation(y1,x1:x2,2))-gapvec;
    else
        gaps = handles.onset(y1:y2,x1:x2);
        gapvec = mean(gaps,1);
        errvec = std(gaps,0,1);
    end
    x = (x1:x2).*handles.xscale;
else
    if get(handles.pointCI,'Value')
        x1 = round((x1+x2)/2);
        gapvec = handles.onset(y1:y2,x1);
        errvec(:,1) = squeeze(handles.onset_variation(y1:y2,x1,1))-gapvec;
        errvec(:,2) = squeeze(handles.onset_variation(y1:y2,x1,2))-gapvec;
        
    else
        [x1,x2,y1,y2] = getarea(handles);
        gaps = handles.onset(y1:y2,x1:x2);
        gapvec = mean(gaps,2);
        errvec = std(gaps,0,2);
    end
    x = (y1:y2).*handles.yscale;
end
figure();
plot(x,gapvec);
if get(handles.errorbarOn,'Value') && sum(sum(abs(errvec)))~=0
    hold on
    if get(handles.pointCI,'Value')
        errorbar(x,gapvec,errvec(:,1),errvec(:,2),'rx');
    else
        errorbar(x,gapvec,errvec,'rx');
    end
    hold off
end
if length(x)>1
    xlim([x(1) x(end)]);
    xlabel(sprintf('%2.2f %s',x(end),handles.xunits),'Fontsize',str2double(get(handles.fontsize,'String')));
end
set(gca,'XTickLabel','','XTick',[],'Fontsize',str2double(get(handles.fontsize,'String')));
ylabel('E_{onset} [eV]','Fontsize',str2double(get(handles.fontsize,'String')));


% Mapping:
function showMap_Callback(hObject, ~, handles)
x = ((1:handles.dim(2))-1).*handles.xscale;
y = ((1:handles.dim(1)) -1).*handles.yscale;
lim1 = str2double(get(handles.statisticsE1,'String'));
lim2 = str2double(get(handles.statisticsE2,'String'));
figure();
imagesc(x,y,handles.onset,[lim1 lim2]);
colormap(get(handles.cmap,'String'));
hcb = colorbar; 
ylabel(hcb,'E_{onset} [eV]')
xlabel(sprintf('%2.2f %s',x(end),handles.xunits),'Fontsize',str2double(get(handles.fontsize,'String')));
ylabel(sprintf('%2.2f %s',y(end),handles.yunits),'Fontsize',str2double(get(handles.fontsize,'String')));
set(gca,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[],'Fontsize',str2double(get(handles.fontsize,'String')));
axis image;
ask = {'Length of scale bar [nm]:'};
answer = inputdlg(ask,'Scale bar?',1,{'10'});
if ~isempty(answer)
    L = str2double(answer{1});
    xlabel('')
    ylabel('')
    set(gca,'XTick',[],'YTick',[],'FontSize',str2double(get(handles.fontsize,'String')))
    N = L/x(2);
    rectangle('Position',[N handles.dim(2)-N N 2],'FaceColor','w','EdgeColor','w')
    str = sprintf('%2.0f nm',N*x(2));
    text(10,83,str,'Color','w','FontSize',str2double(get(handles.fontsize,'String')),'FontWeight','bold')
end

function gapStats_Callback(hObject, ~, handles)
figure();
fig1 = subplot(2,1,1); 
imagesc(handles.onset); 
colormap(get(handles.cmap,'String'));
xlabel('x [pixels]','Fontsize',str2double(get(handles.fontsize,'String')));
ylabel('y [pixels]','Fontsize',str2double(get(handles.fontsize,'String')));
hcb = colorbar; 
ylabel(hcb,'Onset energy [eV]','Fontsize',str2double(get(handles.fontsize,'String')));
h = imrect(fig1,coord2impos(handles));
fcn = makeConstrainToRectFcn('imrect',get(fig1,'XLim'),get(fig1,'YLim'));
setPositionConstraintFcn(h,fcn);
DisplayBandGapWin(getPosition(h),handles)
addNewPositionCallback(h,@(p) DisplayBandGapWin(p,handles));

function DisplayBandGapWin(pos,handles)
x1 = round(pos(1)+0.5);
x2 = round(pos(1)+pos(3)-0.5);
y1 = round(pos(2)+0.5);
y2 = round(pos(4)+pos(2)-0.5);
N = (1+x2-x1)*(1+y2-y1);
g = handles.onset(y1:y2,x1:x2);
R2 = zeros(size(g));
for i=1:size(R2,1)
    for j=1:size(R2,2)
        R2(i,j) = handles.onset_gof{y1+i-1,x1+j-1}.rsquare;
    end
end
[~,p] = vartest(reshape(g,[size(g,1)*size(g,2) 1]),std2(g));
if strcmp(handles.xunits,handles.yunits); unitstr = sprintf('%s^2',handles.xunits);
else; unitstr = sprintf('%s*%s',handles.xunits,handles.yunits);
end

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[50 180 210 20],'fontweight','bold',...
    'String','Onset statistics:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[50 110 130 70],'String',sprintf(...
    'Mean band onset: \nMean standard dev.: \nMaximum onset: \nMinimum onset: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[180 110 80 70],'String',sprintf(...
    ' %1.2f eV\n %1.2f eV\n %1.2f eV\n %1.2f eV',...
    mean2(g),std2(g),max(max(g)),min(min(g)) ));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[50 80 210 20],'fontweight','bold',...
    'String','Region information:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[50 10 90 70],'String',sprintf(...
    'Area: \nPixels: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[140 10 120 70],'String',sprintf(...
    ' %2.3e %s\n [%d:%d;%d:%d] = %d',...
    pos(3)*pos(4)*handles.xscale*handles.yscale,unitstr,x1,x2,y1,y2,N ));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[310 180 210 20],'fontweight','bold',...
    'String','GoF. statistics:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[310 110 130 70],'String',sprintf(...
    'Mean R^2: \nMean standard dev.: \nMaximum R^2: \nMinimum R^2: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[440 110 80 70],'String',sprintf(...
    ' %1.3f eV\n %1.3f eV\n %1.3f eV\n %1.3f eV',...
    mean2(R2),std2(R2),max(max(R2)),min(min(R2))));

uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[310 80 210 20],'fontweight','bold',...
    'String','Distribution:');
uicontrol('style','text','HorizontalAlignment','left',...
    'Position',[310 10 130 70],'String',sprintf(...
    'ANOVA p-value: \nChi-squared p-value: '));
uicontrol('style','text','HorizontalAlignment','right',...
    'Position',[440 10 80 70],'String',sprintf(...
    ' %2.4f\n %1.4f',...
    anova1(g,{},'off'),p));

function savetxt_Callback(hObject, ~, handles)
name = get(handles.origfiletxt,'String');
name = name(1:end-4);
[savename, savepath] = uiputfile('*.txt','Save file as',name);
path = sprintf('%s%s',savepath,savename);
if path ~= 0
    hbar = waitbar(0,'Saving...');
    fileID = fopen(path,'w');
    fprintf(fileID,'Onset energy [eV]:\n');
    for i=1:size(handles.onset,1)
        fprintf(fileID,' %3.5f ',handles.onset(i,:)); fprintf(fileID,'\n');
    end
    fclose(fileID);
    waitbar(1,hbar,'Done'); close(hbar);
end
guidata(hObject,handles)


% Histogram:
function histBins_Callback(hObject, ~, handles)
val = round(str2double(get(hObject,'String')));
val = val*(val>0) + (val<=0);
set(handles.histBins,'String',num2str(val));
guidata(hObject,handles);

function histBins_CreateFcn(hObject, ~, handles)

function gaussianHist_Callback(hObject, ~, handles)

function showHistogram_Callback(hObject, ~, handles)
[x1,x2,y1,y2] = getarea(handles);
onset = handles.onset(y1:y2,x1:x2);

H = reshape(onset,size(onset,1)*size(onset,2),1);
n = str2double(get(handles.histBins,'String'));
figure();
if get(handles.gaussianHist,'Value') == 1
    histfit(H,n)
else
    hist(H,n)
end

title(sprintf('Mean onset: %2.3f eV',mean2(onset)))
xlim([str2double(get(handles.statisticsE1,'String')) str2double(get(handles.statisticsE2,'String'))]);
xlabel('Onset [eV]','Fontsize',str2double(get(handles.fontsize,'String')));
ylabel('Pixel count','Fontsize',str2double(get(handles.fontsize,'String')));
set(gca,'Fontsize',str2double(get(handles.fontsize,'String')));

function getHistogram_Callback(hObject, ~, handles)
figure();
fig1 = subplot(2,1,1);
imagesc(handles.onset);
colormap(get(handles.cmap,'String'));
xlabel('x [pixels]');
ylabel('y [pixels]')
hcb = colorbar;
ylabel(hcb,'Onset energy [eV]')
h = imrect(fig1,coord2impos(handles));
fcn = makeConstrainToRectFcn('imrect',get(fig1,'XLim'),get(fig1,'YLim'));
setPositionConstraintFcn(h,fcn);
DisplayHistogram(getPosition(h),handles)
addNewPositionCallback(h,@(p) DisplayHistogram(p,handles));

function DisplayHistogram(pos,handles)
x1 = round(pos(1)+0.5);
x2 = round(pos(1)+pos(3)-0.5);
y1 = round(pos(2)+0.5);
y2 = round(pos(4)+pos(2)-0.5);
gaps = handles.onset(y1:y2,x1:x2);
subplot(2,1,2);

H = reshape(gaps,size(gaps,1)*size(gaps,2),1);
n = str2double(get(handles.histBins,'String'));
if get(handles.gaussianHist,'Value') == 1
    histfit(H,n)
else
    hist(H,n)
end

title(sprintf('Mean onset: %2.3f eV',mean2(gaps)))
xlabel('Onset [eV]');
ylabel('Pixel count');
xlim([str2double(get(handles.statisticsE1,'String')) str2double(get(handles.statisticsE2,'String'))]);
