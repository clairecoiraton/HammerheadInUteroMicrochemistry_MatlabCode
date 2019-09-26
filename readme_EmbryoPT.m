
% Matlab code used for examining the similarities of transects among embryos of the same mother
% By C. Coiraton May 2017
%All data manipulations and multivariate statistical analyses of were 
%performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
%Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
%College of Marine Science, University of South Florida, St. Petersburg, FL, USA
%http://www.marine.usf.edu/user/djones/matlab/matlab.html 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 INTERPOLATE:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load embryo transect data:
load './2017-03-29/170329.mat' tra*;
load './2017-05-13/170513.mat' tra_59 tra_60;

% Interpolate transects so they are all same length: 
S = {tra_46 tra_47 tra_48 tra_49 tra_50 tra_51 tra_52 tra_53 tra_54 tra_55...
   tra_56 tra_57 tra_58 tra_59 tra_60};    
V = {'tra_46i' 'tra_47i' 'tra_48i' 'tra_49i' 'tra_50i' 'tra_51i' 'tra_52i'...
   'tra_53i' 'tra_54i' 'tra_55i' 'tra_56i' 'tra_57i' 'tra_58i' 'tra_59i' 'tra_60i'};
f_interpolate_PT(S,V);
    
% Clean up:
clear S V tra_46 tra_47 tra_48 tra_49 tra_50 tra_51 tra_52 tra_53 tra_54 tra_55...
   tra_56 tra_57 tra_58 tra_59 tra_60;

% Set up filename & save:
fname = 'Analysis_Embryo_PT';
saver;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   ANALYSIS:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List transects:
tra = (46:60)';

% Create grouping variable (based on transect.xls):
grp = [15 15 15 6 15 6 9 10 6 6 9 10 9 10 6]';

% Show data:
[tra grp]

% Pack interpolated transects into cell array:
S = {tra_46i tra_47i tra_48i tra_49i tra_50i tra_51i tra_52i tra_53i tra_54i...
   tra_55i tra_56i tra_57i tra_58i tra_59i tra_60i};

% Specify target elements (same as in 'readme_Embryo_Mother.m'):
tar = [1 2 3 4 9 10 11 14 17 18 21 22 24];

% Calculate correlation dissimilarity among transects:
d = f_corDis_PT(S,tar);

% Use PCoA to obtain Euclidean coordinates from dissimilarity matrix:
pcoa = f_pcoa(d,0,1,1);
Y    = pcoa.scores;

% Sort the data by group:
[~,key] = sort(grp); % get sort key
Y       = Y(key,:);  % sort Y
grp     = grp(key);  % sort grp

% CAP:
f_capOptimal(Y,'euc',grp,0,1)%centroid
f_capOptimal(Y,'euc',grp,1,1)%spatial median


cap = f_cap(Y,'euc',grp,[],0,1000,1,11,1); % spatial median <- provide optimal m=11
cap = f_cap(Y,'euc',grp,[],1,1000,1,11,1); % centroid <- provide optimal m=11


% Create Plot:

close all;
f_capPlot(cap,{'15' '6' '9' '10'},[],[],[],0.03,'none',0,0,1,0);
figure(1); axis(axis*1.1)
title('Embryo Transects')

% Save plots:
figure(1); f_pdf('Embryo_Transects_CAP')

% Test the significance of the observed classification success rate:
f_chanceClass(grp,1-cap.loo_err.tot,1000,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                RANDOM FOREST:   (works with Matlab R2016b - not 2017                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Random Forests (RF; Breiman 2001) were employed in this study as non-linear 
%counterparts to the CAP classifiers described above. 
%RF were fitted using the same data set and evaluated similarly to estimate reclassification accuracies,
%construct confusion matrices using out-of-bag error rates (OOB-CV) and assess statistical
%significance using PCC. 
%This allowed a comparison of the results and performance of the two distinct methods
%for modeling vertebrae microchemistry data and the selection of the most accurate
%classifier to employ for habitat discrimination (Mercier et al. 2011)

% Generate Random Forest: exploring other method to built better models
% USAGE: model = f_RFclass(X,Y,nTree,mTry,imp,sim,'stnd',verb,X_txt,opt);

rf_GEO = f_RFclass(Y,grp,[],[],0,1,'stnd',1);


% Test the significance of the observed classification success rate:
f_chanceClass(grp,0.7636,1000,1);

%  Create a Random Forest Canonical Discriminant Analysis plot:
f_RFvis(rf_GEO,0.95,1000,0.05);

% Save as PDF:
f_pdf('rf_PT');
f_pdf('rf_PT_vectors');

