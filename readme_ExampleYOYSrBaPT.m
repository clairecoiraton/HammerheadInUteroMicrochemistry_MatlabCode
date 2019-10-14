% Example for finding peaks or valleys in the in utero profile transect data of
% the young-of-the-year captured in Salina Cruz
% By C. Coiraton May 2017

% All data manipulations and multivariate statistical analyses of were 
% performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
% Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
% College of Marine Science, University of South Florida, St. Petersburg, FL, USA
% https://www.marine.usf.edu/research/matlab-resources/fathom-toolbox-for-matlab/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Habitats use changes of mothers during gestation: Sr:Ba 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 91 - SC-LEW-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-05-14/170514.mat tra_91;
X = tra_91;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 92 - SC-LEW-19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-05-14/170514.mat tra_92;
X = tra_92;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TRANSECT 93 - SC-LEW-45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data:
load ./2017-05-14/170514.mat tra_93;
X = tra_93;

% Find valleys in Sr:Ba signals:
[pks,loc] = f_peaks_PT(X,{'Sr88' 'Ba137'},1,0,11,1,3,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
