
% Matlab code used for examining the similarities of transects among embryos of the same mother
% By C. Coiraton May 2017
%All data manipulations and multivariate statistical analyses of were 
%performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
%Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
%College of Marine Science, University of South Florida, St. Petersburg, FL, USA
%http://www.marine.usf.edu/user/djones/matlab/matlab.html 

% List of embryo data from transect.xls:
% 
% date	         loc	num	vert	PT	slide  Embryo
% '29-Mar-2017'	PM	  15	   1	   46	7		 4
% '29-Mar-2017'	PM	  15	   2	   47	7		 2
% '29-Mar-2017'	PM	  15	   3		48	7		 3
% '29-Mar-2017'	PM	  6	   4		49	7		 1
% '29-Mar-2017'	PM	  15	   5		50	7		 1
% '29-Mar-2017'	PM	  6	   6		51	7		 2
% '29-Mar-2017'	PM	  9	   7		52	7		 2
% '29-Mar-2017'	PM	  10	   8		53	7		 1
% '29-Mar-2017'	PM	  6	   9		54	7		 4
% '29-Mar-2017'	PM	  6	   10		55	7		 3
% '29-Mar-2017'	PM	  9	   11		56	7		 3
% '29-Mar-2017'	PM	  10	   12		57	7		 2
% '29-Mar-2017'	PM	  9	   13		58	7		 1
% 13-May-2017'	   PM	  10	   5		59	18		 2
% 13-May-2017'	   PM	  6	   6		60	18		 4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  VARIABLES:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Elements:
% 1   'Li7'
% 2   'Na23'
% 3   'Mg24'
% 4   'P31'
% 5   'Ca43'
% 6   'Sc45'
% 7   'V51'
% 8   'Cr53'
% 9   'Mn55'
% 10  'Fe57'
% 11  'Co59'
% 12  'Ni60'
% 13  'Cu63'
% 14  'Zn64'
% 15  'Cu65'
% 16  'Ge72'
% 17  'Rb85'
% 18  'Sr88'
% 19  'Y89'
% 20  'Cd114'
% 21  'Sn118'
% 22  'Ba137'
% 23  'Au197'
% 24  'Pb208'
% 25  'Th232'
% 26  'U238'


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
f_capOptimal(Y,'euc',grp,0,1)
f_capOptimal(Y,'euc',grp,1,1)


cap = f_cap(Y,'euc',grp,[],0,10000,1,11,1); % <- provide optimal m=11
cap = f_cap(Y,'euc',grp,[],1,10000,1,11,1); % <- provide optimal m=11


% Create Plot:

close all;
%f_capPlot(cap,{'15' '6' '9' '10'},[],[],[],0.03,'none',0,0,0,0);
f_capPlot(cap,{'15' '6' '9' '10'},[],[],[],0.03,'none',0,0,1,0);
figure(1); axis(axis*1.1)
title('Embryo Transects')

% Save plots:
figure(1); f_pdf('Embryo_Transects_CAP')

% Test the significance of the observed classification success rate:
f_chanceClass(grp,1-cap.loo_err.tot,10000,1);

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
f_chanceClass(grp,0.7538,1000,1);


%  Create a Random Forest Canonical Discriminant Analysis plot:
f_RFvis(rf_GEO,0.95,1000,0.05);

% Save as PDF:
f_pdf('rf_PT');
f_pdf('rf_PT_vectors');

