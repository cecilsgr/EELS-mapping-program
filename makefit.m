function [ fitresult ] = makefit( E,I,fitinput )
%MAKEFIT Fit the band onset of a spectrum by the set method
%   Fit the spectrum with energy vector 'E' and intensity vector 
%   'I' using fitting method of number 'fitinput.method_number', 
%   after checking the input parameters of 'fitinput'.
%
%method_number:
%1: Fixed endpoint
%	1.5: Fixed endpoint with dynamic background subtraction
%2: Sliding interval
%	2.5: Sliding interval with dynamic background subtraction
%Other: Regular least squares
%
%fitrange_checking: 1/0, check and display warning for incorrect input parameters

% Enable parameter checking, if not disabled (fitinput.fitrange_checking=0)
if ~isfield(fitinput,'fitrange_checking')   % If field does not exist
    fitinput.fitrange_checking = 1;         % Enable parameter checking
    w=1;
    disp('Fit parameter checking enabled');
else
    w = fitinput.fitrange_checking;
end

% Check fields if enabled:
if w; [E,I,fitinput] = checkfields( E,I,fitinput ); end

% Check ranges and direct to correct fit method:
switch fitinput.fit_method_number
    case 1
        if w;   checkendpointranges( E,fitinput );     end
        fitresult = fixed_endpoint_fit( E,I,fitinput );
    case 1.5
        if w;   checkendpointranges( E,fitinput ); checkBGranges( E,fitinput);	end
        fitresult = fixed_endpoint_background_fit( E,I,fitinput );
    case 2
        if w;   checkintervalranges( E,fitinput );    end
        fitresult = sliding_interval_fit( E,I,fitinput );
    case 2.5
        if w;   checkintervalranges( E,fitinput ); checkBGranges( E,fitinput);	end
        fitresult = sliding_interval_background_fit( E,I,fitinput );
    otherwise
        fitresult = least_squares_fit( E,I,fitinput );
end

end


%% Fitting functions:

function [ fitresult ] = fixed_endpoint_fit( E,I,fitinput )
%ENDPOINTFIT Band gap fitting using a fixed endpoint.
%	Fitting function from band gap to a fixed endpoint. Band gap is the
%	point of maximum R^2 value. 
%   
%fitinput fields:
%	E:                      Vector with spectrum energy, limited by fit range [eV]
%	I:                      Vector with spectrum intensity [a.u]
%	fit_exponent:           Exponent in the fitting: $$I~c(E-E_{onset})^n$$
%   fit_startvalue_const:   Start value of fit constant c
%   testvalue_start:        Lower limit of test onset [eV]
%   testvalue_stop:         Upper limit of test onset [eV]
%	onset_interval:         Not used -- Sliding interval size [eV]
%   onset_endpoint:         Fixed endpoint value [eV]
%   background_range:       Not used -- Size of dynamic background fit window [eV] 
%   background_distance:	Not used -- Distance ignored between test onset and background fit [eV]
%   fit_tolerance:          Tolerance of the fitting at each onset - see Matlab's 'fit' documentation
%   fit_zerointensitycount:	Maximum number of consecutive test values with zero intensity at onset
%   fit_decreasingR2count:	Maximum number of consecutive test values with decreasing R2
%   fit_R2precision:        Precision of R2 value for onset variation
%	show_result:            1/0, Show resulting fit
% 
%fitoutput fields:
%   onset:                  Fitted onset with smallest R^2 [eV]
%   onset_variation:        Calculated variance of onset [eV]
%   onset_test:             Vector with tested onsets [eV]
%   onset_test_R2:          Vector with R2 from fit at test onsets
%   onset_fit:              Fit result at onset
%   onset_gof:              Goodness of fit at onset
%	E_fitted:               Energy vector of fitted onset [eV]
%	I_fitted:               Vector of fitted line [a.u.]
%   background_fit:         Not used -- Background fit at onset
%   background_gof:         Not used -- Goodness of background fit at onset
%   background_R2:          Not used -- Background fit R2
%	E_subtracted:           Not used -- Energy range of background subtracted data
%	I_subtracted:           Not used -- Intensity of background subtracted data


%% Energy limit indices

i1 = findIndex(E,fitinput.testvalue_start,fitinput.fit_tolerance);          % Low gap index
i2 = findIndex(E,fitinput.testvalue_stop,fitinput.fit_tolerance);           % High gap index
i3 = findIndex(E,fitinput.onset_endpoint,fitinput.fit_tolerance);           % High Elim index
Numtests = i2-i1+1;                                                         % Number of test values


%% Store values

onset_test = zeros(Numtests,1);                                             % Tested values
fit_test = cell(Numtests,1);                                                % Fitting output
gof_test = cell(Numtests,1);                                                % Goodness-of-fit
onset_test_R2 = zeros(Numtests+1,1);                                        % R^2 values + extra zero


%% Make fit

fitfunction = fittype( @(c,x) c.*real((x).^(fitinput.fit_exponent)) );      % Fitting function
n_zerointensity = 0;                                                        % Counter for zero intensity
n_decreasingR2 = 0;                                                         % Counter for decreasing R2
i = Numtests;                                                               % Start at highest allowed value

% Iterate from highest test value, as long as counters are below limit: 
while ( i > 0 ) && ( n_zerointensity < fitinput.fit_zerointensitycount ) && ( n_decreasingR2 < fitinput.fit_decreasingR2count )
    % Current test onset:
    onset_test(i) = E(i1+i-1);
    
    % Make fit:
    [fit_test{i},gof_test{i}] = fit(...
        E(i1+i-1:i3)-onset_test(i),I(i1+i-1:i3),...
        fitfunction,'StartPoint',fitinput.fit_startvalue_const,...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    
    % Make valid R2 number:
    onset_test_R2(i) = getR2( gof_test{i}.rsquare );
    
    % Update counters:
    [ n_zerointensity,n_decreasingR2 ] = updatecounts( I(i1+i-1),n_zerointensity,onset_test_R2(i),onset_test_R2(i+1),n_decreasingR2 );
    i = i-1;
end
onset_test_R2(end) = []; % Delete the extra zero in R2


%% Make output

% Find onset from tests:
[~,N2] = max(onset_test_R2);                                                % Index of onset
onset = onset_test(N2);                                                     % Onset value
f = fit_test{N2};                                                           % Fit at onset

% Remove unused elements due to counters:
onset_test(1:i+1) = [];
onset_test_R2(1:i+1) = [];

% Make output fields:
fitresult.onset = onset;                                                    % Fitted onset
fitresult.onset_variation = findvar(onset_test,onset_test_R2,fitinput.fit_R2precision);   % Onset variation
fitresult.onset_fit = f;                                                    % Fit at onset
fitresult.onset_gof = gof_test{N2};                                         % Goodness-of-fit at onset
fitresult.onset_test = onset_test;                                          % Test onsets
fitresult.onset_test_R2 = onset_test_R2;                                    % R^2 at each test onset
fitresult.E_fitted = E(i1+N2-1:i3);                                         % Energy of fit
fitresult.I_fitted = fitfunction(f.c,fitresult.E_fitted-onset);             % Intensity of fit
fitresult.background_fit = 0;                                               % Not in use
fitresult.background_gof = 0;                                               % Not in use
fitresult.background_R2 = 0;                                                % Not in use
fitresult.E_subtracted = E;                                                 % Not in use
fitresult.I_subtracted = I;                                                 % Not in use


%% Show figure

if fitinput.show_result
    plotfig1( E,I,fitresult,fitinput );
end

end


function [ fitresult ] = fixed_endpoint_background_fit( E,I,fitinput )
%ENDPOINTFITBG Band gap fitting and background subtraction using a fixed
%endpoint.
%	Fitting function from band gap to a fixed endpoint, with background
%	subtraction by a power law function. Band gap is the point of maximum
%	R^2 value. 
%   
%fitinput fields:
%	E:                      Vector with spectrum energy, limited by fit range [eV]
%	I:                      Vector with spectrum intensity [a.u]
%	fit_exponent:           Exponent in the fitting: $$I~c(E-E_{onset})^n$$
%   fit_startvalue_const:   Start value of fit constant c
%   testvalue_start:        Lower limit of test onset [eV]
%   testvalue_stop:         Upper limit of test onset [eV]
%	onset_interval:         Not used -- Sliding interval size [eV]
%   onset_endpoint:         Fixed endpoint value [eV]
%   background_range:       Size of dynamic background fit window [eV] 
%   background_distance:	Distance ignored between test onset and background fit [eV]
%   fit_tolerance:          Tolerance of the fitting at each onset - see Matlab's 'fit' documentation
%   fit_zerointensitycount:	Maximum number of consecutive test values with zero intensity at onset
%   fit_decreasingR2count:	Maximum number of consecutive test values with decreasing R2
%   fit_R2precision:        Precision of R2 value for onset variation
%	show_result:            1/0, Show resulting fit
% 
%fitoutput fields:
%   onset:                  Fitted onset with smallest R^2 [eV]
%   onset_variation:        Calculated variance of onset [eV]
%   onset_test:             Vector with tested onsets [eV]
%   onset_test_R2:          Vector with R2 from fit at test onsets
%   onset_fit:              Fit result at onset
%   onset_gof:              Goodness of fit at onset
%	E_fitted:               Energy vector of fitted onset [eV]
%	I_fitted:               Vector of fitted line [a.u.]
%   background_fit:         Background fit at onset
%   background_gof:         Goodness of background fit at onset
%   background_R2:          Background fit R2
%	E_subtracted:           Energy range of background subtracted data
%	I_subtracted:           Intensity of background subtracted data


%% Energy limit indices

dE = E(2)-E(1);
i1 = findIndex(E,fitinput.testvalue_start,fitinput.fit_tolerance);          % Low gap index
i2 = findIndex(E,fitinput.testvalue_stop,fitinput.fit_tolerance);           % High gap index
i3 = findIndex(E,fitinput.onset_endpoint,fitinput.fit_tolerance);           % High Elim index
Numtests = i2-i1+1;                                                         % Number of test values
iBG = i1-round(fitinput.background_range/dE);                               % Background fit index
xBG = round(fitinput.background_distance/dE);                               % Background stop number


%% Store values

onset_test = zeros(Numtests,1);                                             % Tested values
fit_test = cell(Numtests,1);                                                % Fitting output
gof_test = cell(Numtests,1);                                                % Goodness-of-fit
onset_test_R2 = zeros(Numtests+1,1);                                        % R^2 values + extra zero
fit_background = cell(Numtests,1);                                          % Background output
gof_background = cell(Numtests,1);                                          % Background gof
R2_background = zeros(Numtests,1);                                          % R^2 for background


%% Make fit

fitfunction_background = fittype( @(a,r,x) a.*x.^(-r) );                    % Background function
fitfunction = fittype( @(c,x) c.*real((x).^(fitinput.fit_exponent)) );      % Band gap onset function
n_zerointensity = 0;                                                        % Counter for zero intensity
n_decreasingR2 = 0;                                                         % Counter for decreasing R2
i = Numtests;                                                               % Start at highest allowed value

% Iterate from highest test value, as long as counters are below limit: 
while ( i > 0 ) && ( n_zerointensity < fitinput.fit_zerointensitycount ) && ( n_decreasingR2 < fitinput.fit_decreasingR2count )
    % Current test onset:
    onset_test(i) = E(i1+i-1);
    
    % Remove background:
    [fit_background{i},gof_background{i}] = fit(...
        E(iBG+i-1:i1-xBG+i-1),I(iBG+i-1:i1-xBG+i-1),...
        fitfunction_background,'StartPoint',[1 1],...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    R2_background(i) = gof_background{i}.rsquare;
    I2 = I-fitfunction_background(fit_background{i}.a,fit_background{i}.r,E);              % Remove background
    
    % Fit onset:
    [fit_test{i},gof_test{i}] = fit(...
        E(i1+i-1:i3)-onset_test(i),I2(i1+i-1:i3),...
        fitfunction,'StartPoint',fitinput.fit_startvalue_const,...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    
    % Make valid R2 evaluator:
    onset_test_R2(i) = getR2( gof_test{i}.rsquare );
    
    % Update counters:
    [ n_zerointensity,n_decreasingR2 ] = updatecounts( I(i1+i-1),n_zerointensity,onset_test_R2(i),onset_test_R2(i+1),n_decreasingR2 );
    i = i-1;
end
onset_test_R2(end) = []; % Delete the extra zero in R2


%% Make output

% Find onset from tests:
[~,N2] = max(onset_test_R2);                                                % Index of onset
onset = onset_test(N2);                                                     % Onset value
f = fit_test{N2};                                                           % Fit at onset
background_fit = fit_background{N2};                                        % Background fit
background_gof = gof_background{N2};                                        % Background gof

% Remove unused elements due to counters:
onset_test(1:i+1) = [];
onset_test_R2(1:i+1) = [];
R2_background(1:i+1) = [];

% Make output fields:
fitresult.onset = onset;                                                    % Fitted onset
fitresult.onset_variation = findvar(onset_test,onset_test_R2,fitinput.fit_R2precision );   % Onset variation
fitresult.onset_fit = f;                                                    % Fit at onset
fitresult.onset_gof = gof_test{N2};                                         % Goodness-of-fit at onset
fitresult.onset_test = onset_test;                                          % Test onsets
fitresult.onset_test_R2 = onset_test_R2;                                    % R^2 at each test onset
fitresult.E_fitted = E(i1+N2-1:i3);                                         % Energy of fit
fitresult.I_fitted = fitfunction(f.c,fitresult.E_fitted-onset);          	% Intensity of fit
fitresult.background_fit = background_fit;                                  % Background fit
fitresult.background_gof = background_gof;                                  % Background goodness-of-fit
fitresult.background_R2 = R2_background;                                    % Background R^2 at each testgap
fitresult.E_subtracted = E(iBG+N2-1:end);                                   % Inelastic energy
fitresult.I_subtracted = I(iBG+N2-1:end)-fitfunction_background(background_fit.a,background_fit.r,fitresult.E_subtracted); % Inelastic intensity


%% Show figure

if fitinput.show_result
    plotfig2( E,I,fitresult,fitinput,fitfunction_background );
end


end


function [ fitresult ] = sliding_interval_fit( E,I,fitinput )
%SLIDEFIT Band gap fitting using a sliding interval.
%	Fitting function with a sliding fitting interval. Band gap is the point
%	of maximum R^2 value. 
%   
%fitinput fields:
%	E:                      Vector with spectrum energy, limited by fit range [eV]
%	I:                      Vector with spectrum intensity [a.u]
%	fit_exponent:           Exponent in the fitting: $$I~c(E-E_{onset})^n$$
%   fit_startvalue_const:   Start value of fit constant c
%   testvalue_start:        Lower limit of test onset [eV]
%   testvalue_stop:         Upper limit of test onset [eV]
%	onset_interval:         Sliding interval size [eV]
%   onset_endpoint:         Not used -- Fixed endpoint value [eV]
%   background_range:       Not used -- Size of dynamic background fit window [eV] 
%   background_distance:	Not used -- Distance ignored between test onset and background fit [eV]
%   fit_tolerance:          Tolerance of the fitting at each onset - see Matlab's 'fit' documentation
%   fit_zerointensitycount:	Maximum number of consecutive test values with zero intensity at onset
%   fit_decreasingR2count:	Maximum number of consecutive test values with decreasing R2
%   fit_R2precision:        Precision of R2 value for onset variation
%	show_result:            1/0, Show resulting fit
% 
%fitoutput fields:
%   onset:                  Fitted onset with smallest R^2 [eV]
%   onset_variation:        Calculated variance of onset [eV]
%   onset_test:             Vector with tested onsets [eV]
%   onset_test_R2:          Vector with R2 from fit at test onsets
%   onset_fit:              Fit result at onset
%   onset_gof:              Goodness of fit at onset
%	E_fitted:               Energy vector of fitted onset [eV]
%	I_fitted:               Vector of fitted line [a.u.]
%   background_fit:         Not used -- Background fit at onset
%   background_gof:         Not used -- Goodness of background fit at onset
%   background_R2:          Not used -- Background fit R2
%	E_subtracted:           Not used -- Energy range of background subtracted data
%	I_subtracted:           Not used -- Intensity of background subtracted data


%% Energy limit indices

dE = E(2)-E(1);
i1 = findIndex(E,fitinput.testvalue_start,fitinput.fit_tolerance);          % Low gap index
i2 = findIndex(E,fitinput.testvalue_stop,fitinput.fit_tolerance);           % High gap index
L = round(fitinput.onset_interval/dE);                                      % High Elim index
Numtests = i2-i1+1;                                                         % Number of test values


%% Store values

onset_test = zeros(Numtests,1);                                             % Tested onset values
fit_test = cell(Numtests,1);                                               	% Fitting outputs
gof_test = cell(Numtests,1);                                                % Goodness-of-fits
onset_test_R2 = zeros(Numtests+1,1);                                       	% R^2 values, + extra zero


%% Make fit

fitfunction = fittype( @(c,x) c.*real((x).^(fitinput.fit_exponent)) );      % Fitting function
n_zerointensity = 0;                                                        % Counter for zero intensity
n_decreasingR2 = 0;                                                         % Counter for decreasing R2
i = Numtests;                                                               % Start at highest allowed value

% Iterate from highest test value, as long as counters are below limit: 
while ( i > 0 ) && ( n_zerointensity < fitinput.fit_zerointensitycount ) && ( n_decreasingR2 < fitinput.fit_decreasingR2count )
    % Current test onset:
    onset_test(i) = E(i1+i-1);
    
    % Make fit:
    [fit_test{i},gof_test{i}] = fit(...
        E(i1+i-1:i1+L+i-1)-onset_test(i),I(i1+i-1:i1+L+i-1),...
        fitfunction,'StartPoint',fitinput.fit_startvalue_const,...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    
    % Make valid R2 number:
    onset_test_R2(i) = getR2( gof_test{i}.rsquare );
    
    % Update counters:
    [ n_zerointensity, n_decreasingR2 ] = updatecounts( I(i1+i-1),n_zerointensity,onset_test_R2(i),onset_test_R2(i+1),n_decreasingR2 );
    i = i-1;
end
onset_test_R2(end) = []; % Delete the extra zero in R2


%% Make output

% Find onset from tests:
[~,N2] = max(onset_test_R2);                                                % Index of onset
onset = onset_test(N2);                                                     % Onset value
f = fit_test{N2};                                                           % Fit at onset

% Remove unused elements due to counters:
onset_test(1:i+1) = [];
onset_test_R2(1:i+1) = [];

% Make output fields:
fitresult.onset = onset;                                                    % Fitted onset
fitresult.onset_variation = findvar(onset_test,onset_test_R2,fitinput.fit_R2precision);   % Onset variation
fitresult.onset_fit = f;                                                    % Fit at onset
fitresult.onset_gof = gof_test{N2};                                         % Goodness-of-fit at onset
fitresult.onset_test = onset_test;                                          % Test onsets
fitresult.onset_test_R2 = onset_test_R2;                                    % R^2 at each test onset
fitresult.E_fitted = E(i1+N2-1:i1+N2+L-1);                                  % Energy of fit
fitresult.I_fitted = fitfunction(f.c,fitresult.E_fitted-onset);             % Intensity of fit
fitresult.background_fit = 0;                                               % Not in use
fitresult.background_gof = 0;                                               % Not in use
fitresult.background_R2 = 0;                                                % Not in use
fitresult.E_subtracted = E;                                                 % Not in use
fitresult.I_subtracted = I;                                                 % Not in use


%% Show figure

if fitinput.show_result
    plotfig1( E,I,fitresult,fitinput );
end


end


function [ fitresult ] = sliding_interval_background_fit( E,I,fitinput )
%SLIDEFIT_BACKGROUND Band gap fitting and background subtraction using a sliding
%interval.
%	Fitting function with a sliding fitting interval, and subtracting the
%	background by a power law function. Band gap is the point of maximum
%	R^2 value. 
%   
%fitinput fields:
%	E:                      Vector with spectrum energy, limited by fit range [eV]
%	I:                      Vector with spectrum intensity [a.u]
%	fit_exponent:           Exponent in the fitting: $$I~c(E-E_{onset})^n$$
%   fit_startvalue_const:   Start value of fit constant c
%   testvalue_start:        Lower limit of test onset [eV]
%   testvalue_stop:         Upper limit of test onset [eV]
%	onset_interval:         Sliding interval size [eV]
%   onset_endpoint:         Not used -- Fixed endpoint value [eV]
%   background_range:       Size of dynamic background fit window [eV] 
%   background_distance:	Distance ignored between test onset and background fit [eV]
%   fit_tolerance:          Tolerance of the fitting at each onset - see Matlab's 'fit' documentation
%   fit_zerointensitycount:	Maximum number of consecutive test values with zero intensity at onset
%   fit_decreasingR2count:	Maximum number of consecutive test values with decreasing R2
%   fit_R2precision:        Precision of R2 value for onset variation
%	show_result:            1/0, Show resulting fit
% 
%fitoutput fields:
%   onset:                  Fitted onset with smallest R^2 [eV]
%   onset_variation:        Calculated variance of onset [eV]
%   onset_test:             Vector with tested onsets [eV]
%   onset_test_R2:          Vector with R2 from fit at test onsets
%   onset_fit:              Fit result at onset
%   onset_gof:              Goodness of fit at onset
%	E_fitted:               Energy vector of fitted onset [eV]
%	I_fitted:               Vector of fitted line [a.u.]
%   background_fit:         Background fit at onset
%   background_gof:         Goodness of background fit at onset
%   background_R2:          Background fit R2
%	E_subtracted:           Energy range of background subtracted data
%	I_subtracted:           Intensity of background subtracted data


%% Energy limit indices

dE = E(2)-E(1);
i1 = findIndex(E,fitinput.testvalue_start,fitinput.fit_tolerance);          % Low gap index
i2 = findIndex(E,fitinput.testvalue_stop,fitinput.fit_tolerance);           % High gap index
L = round(fitinput.onset_interval/dE);                                      % High Elim index
Numtests = i2-i1+1;                                                         % Number of test values
iBG = i1-round(fitinput.background_range/dE);                               % Background fit index
xBG = round(fitinput.background_distance/dE);                               % Background stop number


%% Store values

onset_test = zeros(Numtests,1);                                             % Tested values
fit_test = cell(Numtests,1);                                                % Fitting output
gof_test = cell(Numtests,1);                                                % Goodness-of-fit
onset_test_R2 = zeros(Numtests+1,1);                                        % R^2 values + extra zero
fit_background = cell(Numtests,1);                                          % Background output
gof_background = cell(Numtests,1);                                          % Background gof
R2_background = zeros(Numtests,1);                                          % R^2 for background


%% Make fit

fitfunction_background = fittype( @(a,r,x) a.*x.^(-r) );                    % Background function
fitfunction = fittype( @(c,x) c.*real((x).^(fitinput.fit_exponent)) );      % Band gap onset function
n_zerointensity = 0;                                                        % Counter for zero intensity
n_decreasingR2 = 0;                                                         % Counter for decreasing R2
i = Numtests;                                                               % Start at highest allowed value

% Iterate from highest test value, as long as counters are below limit: 
while ( i > 0 ) && ( n_zerointensity < fitinput.fit_zerointensitycount ) && ( n_decreasingR2 < fitinput.fit_decreasingR2count )
    % Current test onset:
    onset_test(i) = E(i1+i-1);
    
    % Remove background:
    [fit_background{i},gof_background{i}] = fit(...
        E(iBG+i-1:i1-xBG+i-1),I(iBG+i-1:i1-xBG+i-1),...
        fitfunction_background,'StartPoint',[1 1],...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    R2_background(i) = gof_background{i}.rsquare;
    I2 = I-fitfunction_background(fit_background{i}.a,fit_background{i}.r,E);              % Remove background
    
    % Fit onset:
    [fit_test{i},gof_test{i}] = fit(...
        E(i1+i-1:i1+L+i-1)-onset_test(i),I2(i1+i-1:i1+L+i-1),...
        fitfunction,'StartPoint',fitinput.fit_startvalue_const,...
        'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance);
    
    % Make valid R2 number:
    onset_test_R2(i) = getR2( gof_test{i}.rsquare );
    
    % Update counters:
    [ n_zerointensity,n_decreasingR2 ] = updatecounts( I(i1+i-1),n_zerointensity,onset_test_R2(i),onset_test_R2(i+1),n_decreasingR2 );
    i = i-1;
end
onset_test_R2(end) = []; % Delete the extra zero in R2


%% Make output

% Find onset from tests:
[~,N2] = max(onset_test_R2);                                                % Index of onset
onset = onset_test(N2);                                                     % Onset value
f = fit_test{N2};                                                           % Fit at onset
background_fit = fit_background{N2};                                        % Background fit
background_gof = gof_background{N2};                                        % Background gof

% Remove unused elements due to counters:
onset_test(1:i+1) = [];
onset_test_R2(1:i+1) = [];
R2_background(1:i+1) = [];

% Make output fields:
fitresult.onset = onset;                                                    % Fitted onset
fitresult.onset_variation = findvar(onset_test,onset_test_R2,fitinput.fit_R2precision );   % Onset variation
fitresult.onset_fit = f;                                                    % Fit at onset
fitresult.onset_gof = gof_test{N2};                                         % Goodness-of-fit at onset
fitresult.onset_test = onset_test;                                          % Test onsets
fitresult.onset_test_R2 = onset_test_R2;                                    % R^2 at each test onset
fitresult.E_fitted = E(i1+N2-1:i1+N2+L-1);                                  % Energy of fit
fitresult.I_fitted = fitfunction(f.c,fitresult.E_fitted-onset);          	% Intensity of fit
fitresult.background_fit = background_fit;                                  % Background fit
fitresult.background_gof = background_gof;                                  % Background goodness-of-fit
fitresult.background_R2 = R2_background;                                    % Background R^2 at each testgap
fitresult.E_subtracted = E(iBG+N2-1:end);                                   % Inelastic energy
fitresult.I_subtracted = I(iBG+N2-1:end)-fitfunction_background(background_fit.a,background_fit.r,fitresult.E_subtracted); % Inelastic intensity


%% Show figure

if fitinput.show_result
    plotfig2( E,I,fitresult,fitinput,fitfunction_background );
end

end


function [ fitresult ] = least_squares_fit( E,I,fitinput )
%REGULARFIT Band gap fitting using built-in LLS fitting
%	Fitting function within the full region. Band gap is found using a
%	linear least squares method which is built-in in Matlab
%   
%fitinput fields:
%	E:                      Vector with spectrum energy, limited by fit range [eV]
%	I:                      Vector with spectrum intensity [a.u]
%	fit_exponent:           Exponent in the fitting: $$I~c(E-E_{onset})^n$$
%   fit_startvalue_const:   Start value of fit constant c
%   testvalue_start:        Lower limit of test onset [eV]
%   testvalue_stop:         Upper limit of test onset [eV]
%	onset_interval:         Not used -- Sliding interval size [eV]
%   onset_endpoint:         Not used -- Fixed endpoint value [eV]
%   background_range:       Not used -- Size of dynamic background fit window [eV] 
%   background_distance:	Not used -- Distance ignored between test onset and background fit [eV]
%   fit_tolerance:          Tolerance of the fitting at each onset - see Matlab's 'fit' documentation
%   fit_zerointensitycount:	Not used -- Maximum number of consecutive test values with zero intensity at onset
%   fit_decreasingR2count:	Not used -- Maximum number of consecutive test values with decreasing R2
%   fit_R2precision:        Precision of R2 value for onset variation
%	show_result:            1/0, Show resulting fit
% 
%fitoutput fields:
%   onset:                  Fitted onset with smallest R^2 [eV]
%   onset_variation:        Calculated variance of onset [eV]
%   onset_test:             Not used -- Vector with tested onsets [eV]
%   onset_test_R2:          Not used -- Vector with R2 from fit at test onsets
%   onset_fit:              Fit result at onset
%   onset_gof:              Goodness of fit at onset
%	E_fitted:               Energy vector of fitted onset [eV]
%	I_fitted:               Vector of fitted line [a.u.]
%   background_fit:         Not used -- Background fit at onset
%   background_gof:         Not used -- Goodness of background fit at onset
%   background_R2:          Not used -- Background fit R2
%	E_subtracted:           Not used -- Energy range of background subtracted data
%	I_subtracted:           Not used -- Intensity of background subtracted data


%% Get column vectors

if size(E,2)>size(E,1); E = E'; end                                         % Energy column
if size(I,2)>size(I,1); I = I'; end                                         % Intensity column


%% Make fit

fitfunction= fittype( @(c,g,x) (x>g).*c.*real((x-g).^fitinput.fit_exponent) );    % Fitting function

% Make fit:
[f,gof] = fit(E,I,fitfunction,'StartPoint',[fitinput.fit_startvalue_const 1],...
    'TolFun',fitinput.fit_tolerance,'TolX',fitinput.fit_tolerance,...
    'Lower',[0; E(1)],'Upper',[+Inf; E(end)]);

% Confidence interval:
var = confint(f,1-fitinput.fit_R2precision);                                

%% Make output

% Make output fields:
fitresult.onset = f.g;                                                      % Fitted onset
fitresult.onset_variation = [var(1,2) var(2,2)];                            % Onset variation
fitresult.onset_fit = f;                                                    % Fit at onset
fitresult.onset_gof = gof;                                                  % Goodness-of-fit at onset
fitresult.onset_test = f.g;                                                 % Test onsets
fitresult.onset_test_R2 = gof.rsquare;                                      % R^2 at each test onset
fitresult.E_fitted = E;                                                     % Energy of fit
fitresult.I_fitted = fitfunction(f.c,f.g,E);                                % Intensity of fit
fitresult.background_fit = 0;                                               % Not in use
fitresult.background_gof = 0;                                               % Not in use
fitresult.background_R2 = 0;                                                % Not in use
fitresult.E_subtracted = E;                                                 % Not in use
fitresult.I_subtracted = I;                                                 % Not in use


%% Show figure

if fitinput.show_result
    figure();
    plot(E,I,fitresult.E_fitted,fitresult.I_fitted);
    xlabel('Energy loss [eV]');
    ylabel('Intensity [a.u.]');
    xlim([E(1) E(end)]);
    l = legend('Spectrum','Fit');
    set(l,'location','best','box','off');
    title(sprintf('Fit at E_g=%1.5f (%1.5f,%1.5f) eV',f.g,var(1,2),var(2,2)));
end


end


%% Additional functions:

function n = findIndex( R,x,t )
%FINDINDEX Find first index of value x within vector R, within tolerance t

if x < R(1)                 % Check value within range
    n = 1;
    disp('Warning: A value was chosen outside the energy range!')
elseif x > R(end)           % Check value within range
    n = length(R);
    disp('Warning: A value was chosen outside the energy range!')
elseif x == R(1)            % Check first item
    n = 1;
elseif x == R(end)          % Check last item
    n = length(R);
else                        % Within range of vector
    n = find((R-x)>t,1,'first')-1;
end

end

function [E,I,fitinput] = checkfields( E,I,fitinput )
%CHECKFIELDS Check if fitinput contains all necessary fields

% Make sure energy and intensity are column vectors:
if size(E,2)>size(E,1); E = E'; end
if size(I,2)>size(I,1); I = I'; end

% Check if E and I are compatible:
if length(E) ~= length(I)
	error('Energy and intensity of different length!')
end

% Check each fitinput field.
% If field is not included, set a value and give a warning
if ~isfield(fitinput,'fit_exponent');               fitinput.fit_exponent = 0.5;            disp('Using n=0.5');                                end
if ~isfield(fitinput,'testvalue_start');            fitinput.testvalue_start = E(1);        disp('Setting start test onset to E(1)');           end
if ~isfield(fitinput,'testvalue_stop');             fitinput.testvalue_stop = E(end);       disp('Setting stop test onset to E(end)');          end
if floor(fitinput.fit_method_number) == 1
    if ~isfield(fitinput,'onset_endpoint');         fitinput.onset_endpoint = E(end);       disp('Setting fixed endpoint to E(end)');           end
elseif floor(fitinput.fit_method_number) == 2
    if ~isfield(fitinput,'onset_interval');         fitinput.onset_interval = 0.5;          disp('Setting sliding interval to 0.5 eV');         end
end
if ~isfield(fitinput,'fit_startvalue_const');       fitinput.fit_startvalue_const = 1;      disp('Setting constant start value 1.0');           end
if mod(fitinput.fit_method_number,2) > 0
    if ~isfield(fitinput,'background_range');       fitinput.background_range = 0.5;        disp('Setting background fit range to 0.5 eV');     end
    if ~isfield(fitinput,'background_distance');	fitinput.background_distance = 0;       disp('Background fitted adjacent to onset');        end
end
if ~isfield(fitinput,'fit_tolerance');              fitinput.fit_tolerance = 1E-6;          disp('Using default tolerance of Matlab');          end
if ~isfield(fitinput,'fit_zerointensitycount');     fitinput.fit_zerointensitycount = +Inf; disp('No limit of cons. zero int. at test onset');	end
if ~isfield(fitinput,'fit_decreasingR2count');      fitinput.fit_decreasingR2count = +Inf;  disp('No limit of consecutive decreasing R2');      end
if ~isfield(fitinput,'fit_R2precision');            fitinput.fit_R2precision = 0.05;        disp('Using R2 precision of 0.05');                 end
if ~isfield(fitinput,'show_result');                fitinput.show_result = 1;               disp('Showing fit result');                         end

end

function checkintervalranges( E,fitinput )
%CHECKSLIDERANGES Check if fitinput parameters are compatible with the
%sliding interval fit method

if fitinput.testvalue_stop + fitinput.onset_interval - E(end) > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s\n%s','Highest band gap fit higher than energy range',...
        'Solutions: ',...
        '   1) Increase high fit range',...
        '   2) Decrease sliding range (fitinput.interval)',...
        '   3) Decrease upper gap limit (fitinput.fitinput.testvalue_stop)');
end
if E(1) - fitinput.testvalue_start > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s','Lowest band gap fit lower than energy range',...
        'Solutions: ',...
        '   1) Decrease low fit range',...
        '   2) Increase lower gap limit (fitinput.testvalue_start)');
end
end

function checkendpointranges( E,fitinput )
%CHECKENTPTRANGES Check if fitinput parameters are compatible with the
%fixed endpoint fit method

if fitinput.onset_endpoint - E(end) > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s','Endpoint not in energy vector!',...
        'Solutions: ',...
        '   1) Increase high fit range',...
        '   2) Decrease endpoint value (fitinput.onset_endpoint)');
end
if fitinput.testvalue_stop - E(end) > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s','Highest band gap fit higher than energy range',...
        'Solutions: ',...
        '   1) Increase high fit range',...
        '   2) Decrease upper gap limit (fitinput.testvalue_stop)');
end
if fitinput.testvalue_stop + (E(3)-E(1)) - fitinput.onset_endpoint > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s','Too few points fitted at high gap values',...
        'Solutions: ',...
        '   1) Increase endpoint (fitinput.onset_endpoint)',...
        '   2) Decrease upper gap limit (fitinput.testvalue_stop)');
end
if E(1) - fitinput.testvalue_start > fitinput.fit_tolerance
    error('%s\n%s\n%s\n%s','Lowest band gap fit lower than energy range',...
        'Solutions: ',...
        '   1) Decrease low fit range',...
        '   2) Increase lower gap limit (fitinput.testvalue_start)');
end
end

function checkBGranges( E,fitinput)

if fitinput.background_range < fitinput.background_distance
    error('%s\n%s\n%s\n%s','Background fit larger than distance',...
        'Solutions: ',...
        '   1) Increase background fit range (fitinput.background_range)',...
        '   2) Decrease distance without background subtraction (fitinput.background_distance)');
end

if E(1) > fitinput.testvalue_start - fitinput.background_range
    error('%s\n%s\n%s\n%s\n%s','Lowest band gap fit background range outside energy range',...
        'Solutions: ',...
        '   1) Decrease low fit range',...
        '   2) Increase background fit range (fitinput.background_range)',...
        '   3) Increase lower gap limit (fitinput.testvalue_start)');
end
end

function r2 = getR2( r2 )
%getR2 return R2 as a valid number
if isfinite(r2)         % R^2 is number
    r2 = r2*(r2>0);     % R^2 is real
else
    r2 = 0;             % ... or zero
end
end

function [ zerointensitycount,decreasingR2count ] = updatecounts( intensity,zerointensitycount,R2,previousR2,decreasingR2count )
%UPDATECOUNTS Update counters
if intensity<=0
    zerointensitycount = zerointensitycount+1;
else
    zerointensitycount = 0;
end
if R2<previousR2
    decreasingR2count = decreasingR2count+1;
else
    decreasingR2count = 0;
end

end

function vars = findvar( testonset,R2,R2diff )
% FINDVAR Find closest bounds to intersection between R2 and R2precision
[Rm,N2] = max(R2);                                          % Max R^2
y = Rm-R2diff;                                              % Intersection
j = find(diff(R2 <= y));                                    % Intersection pts

if isempty(j)
    vars = [testonset(N2) testonset(N2)];
elseif size(j) == 1
    if j(1) < N2
        vars = [testonset(j(1)) testonset(N2)];
    else
        vars = [testonset(N2) testonset(j(1))];
    end
else
    if j(1)>1 && j(end)<length(testonset)
        bounds = testonset(j)+(y-R2(j)).*(testonset(j+1)-testonset(j))./(R2(j+1)-R2(j));
        if size(bounds,1) > 2
            bounds = [bounds(find(bounds<testonset(N2),1,'last')) bounds(find(bounds>testonset(N2),1))];
        end
        vars = [min(bounds) max(bounds)];
    else
        if j(1) == 1
            vars(1) = testonset(1);
        else
            vars(1) = testonset(j);
        end
        if j(end) == length(testonset)
            vars(2) = testonset(end);
        else
            vars(2) = testonset(j(end));
        end
    end
end
% if size(testonset) < j(end)+1
%     bounds = [0 0];
%     disp('a')
% elseif j(1) < 1
%     bounds = [0 0];
%     disp('b')
% else
%     bounds = testonset(j)+(y-R2(j)).*(testonset(j+1)-testonset(j))./(R2(j+1)-R2(j));
% end

%if size(bounds,1) == 0; bounds = [0 0]; end
% if size(bounds,1) > 2
%     bounds = [bounds(find(bounds<testonset(N2),1,'last')) bounds(find(bounds>testonset(N2),1))];
% end
% vars = [min(bounds) max(bounds)];
end

function plotfig1( E,I,fitresult,fitinput )
%PLOTFIG1 Plot onset fit

figure();

subplot(1,2,1)
plot(E,I,fitresult.E_fitted,fitresult.I_fitted);
xlabel('Energy loss [eV]');
ylabel('Intensity [a.u.]');
xlim([fitresult.E_fitted(1)-0-1 fitresult.E_fitted(end)+0.1]);
limsy=get(gca,'YLim');
ylim([0 limsy(2)])
l = legend('Spectrum','Fit');
set(l,'location','best','box','off');
title(sprintf('Fit at E_g=%1.2f eV',fitresult.onset));

subplot(1,2,2)
Rm = max(fitresult.onset_test_R2);
y = Rm-fitinput.fit_R2precision;

plot(fitresult.onset_test,fitresult.onset_test_R2,fitresult.onset,Rm,'k*')
if fitresult.onset_variation(1)>0 && fitresult.onset_variation(2)>0
    hold on
    plot([fitresult.onset_test(1) fitresult.onset_test(end)],[y y],'k--',...
        fitresult.onset_variation,[y,y],'k*',...
        [fitresult.onset_variation(1) fitresult.onset_variation(1)],[y 0],'k--',...
        [fitresult.onset_variation(2) fitresult.onset_variation(2)],[y 0],'k--');
    title(sprintf('Onset: %1.2f (%1.2f, %1.2f) eV',fitresult.onset,fitresult.onset_variation));
else
    title(sprintf('Onset: %1.2f eV',fitresult.onset));
end
xlabel('E_g [eV]');
ylabel('R^2 [a.u.]');
xlim([fitresult.onset_test(1) fitresult.onset_test(end)]);

end

function plotfig2( E,I,fitresult,fitinput,fitfunction_background )
%PLOTFIG2 Plot background fit and onset fit

ind = findIndex(fitresult.E_subtracted,fitresult.onset,0); 
xBG = round(fitinput.background_distance/E(2)-E(1));

figure();

subplot(2,2,1)
hold on
a = area(fitresult.E_subtracted(1:(ind-xBG)),...
   fitfunction_background(fitresult.background_fit.a,fitresult.background_fit.r,fitresult.E_subtracted(1:(ind-xBG))));
p = plot(E,I,fitresult.E_subtracted,...
    fitfunction_background(fitresult.background_fit.a,fitresult.background_fit.r,fitresult.E_subtracted),...
    fitresult.E_subtracted,fitresult.I_subtracted);
hold off
xlabel('Energy loss [eV]');
ylabel('Intensity [a.u.]');
xlim([fitresult.E_subtracted(1) fitresult.E_fitted(end)+0.1]);
limsy=get(gca,'YLim');
ylim([0 limsy(2)])
l = legend(p,'Spectrum','Background','Inelastic');
set(l,'location','best','box','off');
set(a,'FaceColor',[1 0.5 0.4],'EdgeColor','none');
title(sprintf('Background subtraction\n at E_{onset}=%1.2f eV',fitresult.onset));

subplot(2,2,3)
plot(fitresult.onset_test,fitresult.background_R2,fitresult.onset,fitresult.background_gof.rsquare,'k*')
xlabel('E_{onset} [eV]');
ylabel('Background R^2 [a.u.]');
xlim([fitresult.onset_test(1) fitresult.onset_test(end)]);
title('Background fit R^2');

subplot(2,2,2)
plot(fitresult.E_subtracted,fitresult.I_subtracted,fitresult.E_fitted,fitresult.I_fitted);
xlabel('Energy loss [eV]');
ylabel('Intensity [a.u.]');
xlim([fitresult.E_fitted(1)-fitinput.background_range-0.1 fitresult.E_fitted(end)+0.1]);
limsy=get(gca,'YLim');
ylim([0 limsy(2)])
l = legend('Spectrum','Fit');
set(l,'location','best','box','off');
title(sprintf('Fit at E_{onset}=%1.2f eV',fitresult.onset));

subplot(2,2,4)
Rm = max(fitresult.onset_test_R2);
y = Rm-fitinput.fit_R2precision;
plot(fitresult.onset_test,fitresult.onset_test_R2,fitresult.onset,Rm,'k*')
if fitresult.onset_variation(1)>0 && fitresult.onset_variation(2)>0
    hold on
    plot([fitresult.onset_test(1) fitresult.onset_test(end)],[y y],'k--',...
        fitresult.onset_variation,[y,y],'k*',...
        [fitresult.onset_variation(1) fitresult.onset_variation(1)],[y 0],'k--',...
        [fitresult.onset_variation(2) fitresult.onset_variation(2)],[y 0],'k--');
    title(sprintf('Onset: %1.2f (%1.2f, %1.2f) eV',fitresult.onset,fitresult.onset_variation));
else
    title(sprintf('Onset: %1.2f eV',fitresult.onset));
end
xlabel('E_{onset} [eV]');
ylabel('R^2 [a.u.]');
xlim([fitresult.onset_test(1) fitresult.onset_test(end)]);

end


