% Example for finding peaks or valleys in the in utero profile transect data of
% embryos
%By C. Coiraton May 2017
%All data manipulations and multivariate statistical analyses of were 
%performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
%Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
%College of Marine Science, University of South Florida, St. Petersburg, FL, USA
%http://www.marine.usf.edu/user/djones/matlab/matlab.html 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Habitats use changes (Sr:Ba) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 46 - PM-LEW-15-E4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-29/170329.mat tra_46;
X = tra_46;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,50,1,10,10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 47 - PM-LEW-15-E2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-29/170329.mat tra_47;
X = tra_47;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,50,1,10,10);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 48 - PM-LEW-15-E3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-03-29/170329.mat tra_48;
X = tra_48;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,50,1,10,10);



