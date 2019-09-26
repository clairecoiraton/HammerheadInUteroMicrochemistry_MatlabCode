% Matlab code used for comparing the elemental signatures deposited at the
% vertebral edge between each pregnant female and her embryos.
%By C. Coiraton May 2017

%All data manipulations and multivariate statistical analyses of were 
%performed using the free download Fathom Toolbox for MatlabTM (Jones DL, 2017)
%Fathom Toolbox for Matlab: software for multivariate ecological and oceanographic data analysis.
%College of Marine Science, University of South Florida, St. Petersburg, FL, USA
%http://www.marine.usf.edu/user/djones/matlab/matlab.html 



% Load data:
load all_HH.mat HH spots;

% Set up filename & save:
fname = 'HH_analysis_Embryo';
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
%   .type   = type of spot sample (FO, BM, GB, ED)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 EDGE SAMPLES:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get index to Embryos:
idxE = (spots.age==-1);

% Get index to spot samples of ED (type_Rx = CO=1, PB=2, BM=3, GB=4, ED=5)
idxED = (spots.type_Rx==5);

% Get index Embryo-ED samples:
idxE_ED = logical(idxE .* idxED);

% Extract data:
Y   = HH.ppm(idxE_ED,:);
LOD = HH.LOD(idxE_ED,:);

% Get percentage of each element below LOD:
bin  = Y<=LOD;
pLOD = (sum(bin)/size(bin,1));
pLOD'

% Get index to target elements (= not consistently below LOD):
tar = find(pLOD<0.1);
tar(ismember(tar,[5])) = []; % remove Ca43 
HH.txt_iso(tar)'                % show target elements

% Confirm target elements are above LOD (should be all 0's):
Y(:,tar)<=LOD(:,tar)

% Get index to MOTHER-ED samples:
idxPM  = spots.loc_Rx==13;
idxNum = ismember(spots.num,[6 9 10 15]);
idxED  = spots.ED_Rx==1;
idxF   = ~(spots.age==-1);
idxMom = logical(idxPM .* idxNum .* idxED .* idxF);
sum(idxMom)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      ANALYSIS - COMBINE EMBRYO & MOTHER:                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get index to Embryo + Mother:
idxEM = logical(idxE_ED + idxMom);

% Extract target elements:
Y   = HH.ppm(idxEM,tar);
grp = spots.num(idxEM);
age = spots.age(idxEM);

% Sort the data by group:
[~,key] = sort(grp); % get sort key
Y       = Y(key,:);  % sort Y
grp     = grp(key);  % sort grp
age     = age(key);  % sort age


% CAP:
f_capOptimal(Y,'euc',grp,1,1);%spatial median

cap = f_cap(Y,'euc',grp,[],1,1000,1,12,1); % <- provide optimal m=12

% Create plots:
close all;
%[~,~,~,crds] = f_capPlot(cap,f_unique(f_num2cell(grp)),[],Y,HH.txt(tar),0.05,'none',0,0,0,0);
[~,~,~,crds] = f_capPlot(cap,f_unique(f_num2cell(grp)),[],Y,HH.txt(tar),0.05,'none',0,0,1,0);
figure(1); axis(axis*1.1); title('Embryo EDGE CAP');
figure(2); title('Embryo EDGE VECTORS');

% Identify mothers:
figure(1);
idxM = (age>0); % index to just mothers
hold on;
plot(crds(idxM,1),crds(idxM,2),'o','MarkerFaceColor','none',...
   'MarkerEdgeColor','k', 'MarkerSize',16)
hold off;

% Save plots:
figure(1); f_pdf('Embryo_ED_CAP')
figure(2); f_pdf('Embryo_ED_VECTORS')

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
f_pdf('rf_Embryos&MomED');
f_pdf('rf_Embryos&MomED_vectors');

