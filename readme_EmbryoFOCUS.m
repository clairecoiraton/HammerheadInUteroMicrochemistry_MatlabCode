
% Matlab code used for comparing the elemental signatures deposited at the
% vertebral focus of the embryos from each litter.
%By C. Coiraton May 2017
%All data manipulations and multivariate statistical analyses of were 
%performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
%Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
%College of Marine Science, University of South Florida, St. Petersburg, FL, USA
%http://www.marine.usf.edu/user/djones/matlab/matlab.html 



% List of embryo data from transect.xls:
% 
% date	         loc	num	vert	PT	slide  Embryo
% '29-Mar-2017'	PM	  15   1	   46	7		 4
% '29-Mar-2017'	PM	  15   2	   47	7		 2
% '29-Mar-2017'	PM	  15 . 3		48	7		 3
% '29-Mar-2017'	PM	  6	   4		49	7		 1
% '29-Mar-2017'	PM	  15   5		50	7		 1
% '29-Mar-2017'	PM	  6	   6		51	7		 2
% '29-Mar-2017'	PM	  9	   7		52	7		 2
% '29-Mar-2017'	PM	  10   8		53	7		 1
% '29-Mar-2017'	PM	  6	   9		54	7		 4
% '29-Mar-2017'	PM	  6	   10		55	7		 3
% '29-Mar-2017'	PM	  9	   11		56	7		 3
% '29-Mar-2017'	PM	  10   12		57	7		 2
% '29-Mar-2017'	PM	  9	   13		58	7		 1
% 13-May-2017'	PM	  10   5		59	18		 2
% 13-May-2017'	PM	  6	   6		60	18		 4
% 
% - there are n = 15 samples
% - includes CO, ED, and PT
% 
% 170329: no FO to match the ED sample for spot 580, must've been broken?
% 170329: no ED to match the FO sample for spot 585, must've been broken? 

% Load data:
load all_HH.mat HH spots;

% Set up filename & save:
fname = 'HH_analysis_Embryo_Core';
saver;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  VARIABLES:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HH = structure with the following fields:
%    .oto     = numeric tag identifying individual otoliths
%    .txt     = cell array of text for each analyte element
%    .iso     = corresponding isotope number
%    .txt_iso = combined element + isotope labels
%    .nRep    = # within-otolith replicate spot samples used to calclate averages
%    .ppm     = concentration of unknown, averaged across replicates                 (ppm)
%    .LOD     = limits of detection, averaged across replicates                      (ppm)
%    .ratio   = molar ratios to internal standard, averaged across replicates (mMole/Mole)
%    .SRM     = name of SRM used for external calibration
%    .adj     = type of adjustment applied to values below LOD             (zero,LOD,none)
%    .spike   = type of spike removal applied to the time series
%    .drift   = method of drift correction
%    .tol     = tolerance for linear interpolation
%    .SRM     = name of SRM used for external calibration
%    .gDate   = cell array of Gregorian date of acquisition
% 
% 
% spots = structure with the following fields:
%   .spot   = unique spot scan identifier (average of 2-3 replicates)
%   .type   = type of spot sample (CO, PB, BM, GB, ED)
%   .loc    = location of collection site
%   .num    = specimen number from this site
%   .yr     = collection year
%   .mo     = collection month
%   .age    = age in years (-1 = embryo)
%   .sex    = F = female, M = male, E = embryo
%   .embryo = embryo number (0 = adult)
%   .slide  = master slide #
%   .vert   = number of vertebrae thin section on master slide
%   .gDate  = date of LA-ICP-MS run
% 
%   .loc_Rx  = location code
%   .type_Rx = (FO=1, PB=2, BM=3, GB=4, ED=5)
%   .GB_Rx   = growth band code (0-14)
%   .ED_Rx   = edge code (0-1)

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
%                                 FOCUS SAMPLES:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get index to Embryos:
idxE = (spots.age==-1);

% Get index to spot samples of FO (type_Rx = (FO=1, PB=2, BM=3, GB=4, ED=5)
idxFO = (spots.type_Rx==1);

% Get index Embryo-CO samples:
idxE_FO = logical(idxE .* idxFO);

% Extract data:
Y   = HH.ppm(idxE_FO,:);
LOD = HH.LOD(idxE_FO,:);

% Get percentage of each element below LOD:
bin  = Y<=LOD;
pLOD = (sum(bin)/size(bin,1));
pLOD'

% Get index to target elements (= not consistently below LOD): Give us a list of
% which elements are the target
tar = find(pLOD<0.1)'
tar(ismember(tar,[5 15]))=[]; % additionally, remove Ca43 (int.std) & Cu65 (Cu63 alr.) [5 15]
HH.txt_iso(tar)'

% Confirm target elements are above LOD (should be all 0's):
Y(:,tar)<=LOD(:,tar)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                ANALYSIS - FOCUS:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract target elements:
Y   = HH.ppm(idxE_FO,tar);
grp = spots.num(idxE_FO);

% Sort the data by group:
[~,key] = sort(grp); % get sort key
Y       = Y(key,:);  % sort Y
grp     = grp(key);  % sort grp

% Create Euclidean distance matrix from standardized data:
yDis = f_dis(f_stnd(Y),'euc');

% CAP:
f_capOptimal(Y,'euc',grp,1,1);%spatial median
f_capOptimal(Y,'euc',grp,0,1);%centroid


cap = f_cap(Y,'euc',grp,[],1,1000,1,9,1);% <- provides optimal m=9 SPATIAL MEDIAN
cap = f_cap(Y,'euc',grp,[],0,1000,1,9,1); % <- provides optimal m=9 CENTROID


% Create plots:
close all; % close previous plots

%f_capPlot(cap,f_unique(f_num2cell(grp)),[],Y,HH.txt(tar),0.05,'none',0,0,0,0);
f_capPlot(cap,f_unique(f_num2cell(grp)),[],Y,HH.txt(tar),0.05,'none',0,0,1,0);
figure(1); axis(axis*1.1); title('Embryo FOCUS CAP');
figure(2); title('Embryo FOCUS VECTORS');

% Save plots:
figure(1); f_pdf('Embryo_FO_CAP')
figure(2); f_pdf('Embryo_FO_VECTORS')

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
f_chanceClass(grp,0.7538,1000,1);


%  Create a Random Forest Canonical Discriminant Analysis plot:
f_RFvis(rf_GEO,0.95,1000,0.05);

% Save as PDF:
f_pdf('rf_EmbryosFO');
f_pdf('rf_EmbryosFO_vectors');

