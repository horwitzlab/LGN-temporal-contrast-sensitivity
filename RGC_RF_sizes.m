% Getting formulas for LGN RF size as a function of retinal eccentricity to
% be used in IsoSampPop.m
% Contents:
% Section 0: Generic constants
% Section 1: Perry et al. 1984 parasol RGC data (anatomy)
% Section 2: Perry et al. 1984 midget RGC data (anatomy)
% Section 3: Croner and Kaplan parasol RGC data (electrophysiology)
% Section 4: Croner and Kaplan midget RGC data (electrophysiology)
% Section 5: Dacey and Petersen midget RGC data (anatomy)
% Section 6: Derrington and Lennie 1984 parvocellular LGN data (electrophysiology)
% Section 7: Derrington and Lennie 1984 magnocellular LGN data (electrophysiology)
% Section 8: Watson 2014 midget RGC data (anatomy + model)
% Section 9: Watanbe and Rodieck 1989 parasol RGC data (anatomy)
% Section 10: Watanbe and Rodieck 1989 midget RGC data (anatomy)
% Section 11: Cone inverse densities

% Note: Temporal retina has fewer cones than nasal retina

%%
% Section 0
% rf_r_mm is the RF eccentricity of the RF center. At some point will have
% to deal with the difference between nasal and temporal eccentricity (r).
MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1./MMPERDEG;

%% Section 1
% Parasol RGC dentritic fields from Perry et al 1984, Figure 6A. 
% Captured with PlotDigitizer.
% Only digitizing data from the central 3 mm = 13.4 deg (there's a gap in
% the data after this and the data appear to get noiser).
% The difference between nasal and temporal RF starts to become apparent
% at eccentricities > ~4 mm (~16°)?

% Croner and Kaplan Figure 13A suggest (with slightly thin data) that
% their physiological measurements of M-cell RF centers are predicted by
% M-cell dendritic field sizes.

perryMdata = xlsread ('RGC_RFsizes_vs_eccentricity','PerryEtAl1984ParasolData'); % x = eccentricity (mm), y = dendritic field diam. (microns)

% Comparing captured data to original. Compare to Figure 6(A).
figure; axes; hold on;
plot(perryMdata(:,1),perryMdata(:,2),'ko');
title('Perry et al. parasol Figure 6A');
xlabel('Eccentricity (mm)');
ylabel('Dendritic field diam. (microns)');

% Comparing to Croner and Kaplan Figure 13a
% Dividing 'y' by 2 to convert from RF diameter to RF radius.
% Dividing 'y' by 1000 to convert  RF diameter from microns to mm.
% RFs are a little too small.
figure; subplot; hold on; 
plot(perryMdata(:,1)*DEGPERMM,perryMdata(:,2)*DEGPERMM/1000/2,'ko');
x = [0 10];
b = regress(perryMdata(:,2)*DEGPERMM/1000/2,[ones(size(perryMdata,1),1) perryMdata(:,1)*DEGPERMM]); % intercept and slope in degrees (DF diam) as a function of degrees (eccentrcity)
plot(x,[[1 x(1)]*b [1 x(2)]*b])
ylabel('RF diameter (deg)');
xlabel('RF eccentricity (deg)')
set(gca,'Ytick',[0.05 .1 .15 .2 .25])
title('Perry et al. parasol for comparison with Croner and Kaplan');

% To compare across studies I will try to plot RF diameter (1 SD in deg) as a
% function of eccentricity (in deg) for every data set.
% Even assuming that the DF is 1 SD of the RF (which seems bizarre) the
% Perry estimates of RF sizes are smaller than Croner and Kaplan's.

figure; axes; hold on;
plot(perryMdata(:,1)*DEGPERMM,perryMdata(:,2)*DEGPERMM/1000,'ko');
x = linspace(0,18,100);
b = regress(perryMdata(:,2)*DEGPERMM/1000,[ones(size(perryMdata,1),1) perryMdata(:,1)*DEGPERMM]); % intercept and slope in degrees (DF diam) as a function of degrees (eccentrcity)
plot(x,[ones(length(x),1) x']*b);
ylabel('RF diameter (deg)');
xlabel('RF eccentricity (deg)')
set(gca,'Ytick',[0:.05:.35])
title('Perry et al. parasol full DF diameter');
set(gca,'Xlim',[0 18],'Ylim',[0 0.35])


%% Section 2 Perry et al. midget data
% Looking at the midget RGC DF diameters
perryPdata = xlsread ('RGC_RFsizes_vs_eccentricity','PerryEtAl1984MidgetData'); % x = eccentricity (mm), y = dendritic field diam. (microns)
figure; axes; hold on;
plot(perryPdata(:,1),perryPdata(:,2),'ko');
% looks like the log of DF diameter increases with eccentricity
set(gca,'Yscale','log');
title('Perry et al. midget Figure 6B');
xlabel('Eccentricity (mm)');
ylabel('Dendritic field diam. (microns)');

% Now converting from x = eccentricity (mm), y = dendritic field diam. (microns)
% to x = eccentricity (deg), y = dendritic field diam. (deg)
figure; axes; hold on;
plot(perryPdata(:,1)*DEGPERMM,perryPdata(:,2)*DEGPERMM/1000,'ko');
b = regress(log10(perryPdata(:,2)*DEGPERMM/1000),[ones(size(perryPdata,1),1) perryPdata(:,1)*DEGPERMM]);
x = linspace(0,18,100);
plot(x,10.^([ones(length(x),1) x']*b));
ylabel('RF diameter (deg)'); % Note this is the full diameter, not 1 SD
xlabel('RF eccentricity (deg)');
title('Perry et al. midget full DF diameter');
set(gca,'Ylim',[0 .2],'Xlim',[0 18])
% This might slightly overestimate the size of RFs at small eccentricities?

%% Section 3
% Croner and Kaplan parasol RGC data
CK_Pdata = xlsread ('RGC_RFsizes_vs_eccentricity','CronerKaplan1994ParasolData'); % x = eccentricity (temporal equivalent ecc. (deg)), y = center radius (deg)
L = CK_Pdata(:,1) < 20; % Only looking at neurons with RFs within the central 20°

% Compare to Croner and Kaplan 1994 Figure 4A open symbols
figure; axes; hold on;
plot(CK_Pdata(L,1),CK_Pdata(L,2),'ko');
b = regress(CK_Pdata(L,2),[ones(sum(L),1) CK_Pdata(L,1)]);
x = linspace(0,18,100);
plot(x,[ones(length(x),1) x']*b);
ylabel('Center radius (deg)'); xlabel('Temporal equivalent eccentricity (deg)');
title('Croner and Kaplan parasol. Compare to Figure 4A open symbols');

% To compare across studies I will try to plot RF diameter (1 SD in deg) as a
% function of eccentricity (in deg) for every data set.
% "*2" to convert radius to diameter (diameter is 2 SDs).
% "/sqrt(2)" to go out only 1 SD instead of 1/e
figure; hold on;
plot(CK_Pdata(L,1),CK_Pdata(L,2)*2/sqrt(2),'ko');
b = regress(CK_Pdata(L,2)*2/sqrt(2),[ones(sum(L),1) CK_Pdata(L,1)]);
x = linspace(0,18,100);
plot(x,[ones(length(x),1) x']*b);
ylabel('Center diameter (1 SD in deg)'); xlabel('Temporal equivalent eccentricity (deg)');
title('Croner and Kaplan parasol. Standard parameterization (1 SD).');
set(gca,'Xlim',[0 18],'Ylim',[0 .35]);
%%
% Section 4: Croner and Kaplan midget data

CK_midget_data = xlsread('RGC_RFsizes_vs_eccentricity','CronerKaplan1994MidgetData'); % x = eccentricity (temporal equivalent ecc. (deg)), y = center radius (deg)
L = CK_midget_data(:,1) < 20; % Only looking at neurons with RFs within the central 20°

% Compare to Croner and Kaplan 1994 Figure 4A open symbols
figure; axes; hold on;
plot(CK_midget_data(L,1),CK_midget_data(L,2),'ko');
b = regress(log10(CK_midget_data(L,2)),[ones(sum(L),1) CK_midget_data(L,1)]);
x = linspace(0,18,100);
plot(x,10.^([ones(length(x),1) x']*b));
ylabel('Center radius (deg)'); xlabel('Temporal equivalent eccentricity (deg)');
title('Croner and Kaplan midget for comparison with Figure 4A');

% To compare across studies I will try to plot RF diameter (1 SD in deg) as a
% function of eccentricity (in deg) for every data set.
% "*2" to convert radius to diameter (diameter is 2 SDs)
% "/sqrt(2)" to go out only 1 SD instead of 1/e
figure; hold on;
plot(CK_midget_data(L,1),CK_midget_data(L,2)*2/sqrt(2),'ko');
b = regress(log10(CK_midget_data(L,2)*2/sqrt(2)),[ones(sum(L),1) CK_midget_data(L,1)]);
x = linspace(0,18,100);
plot(x,10.^([ones(length(x),1) x']*b));
ylabel('Center diameter (1 SD in deg)'); xlabel('Temporal equivalent eccentricity (deg)');
title('Croner and Kaplan midget 1 SD');
set(gca,'Ylim',[0 .2],'Xlim',[0 18])

%%
% Section 5: Dacey and Petersen 1992 (midgets only)
% Legend to Figure 2 gives an equation for the midget dendritic field size
% (in microns per mm from fovea for temporal, upper, and lower retina) 
% Dacey's equation (y = 8.64*x.^1.04 from Figure 2A) doesn't fit his data
% well over the range that I care about. 

DP_midget_data = xlsread ('RGC_RFsizes_vs_eccentricity','DaceyPetersen1992MidgetData'); % x = eccentricity (deg), y = dentritic field diam. (arc min)
L = DP_midget_data(:,1) < 20;
figure; axes; hold on;
plot(DP_midget_data(L,1),DP_midget_data(L,2),'k^');
% Compare to Figure 2B

% /60 to get from arc mins to degrees
b = regress(log10(DP_midget_data(L,2)/60),[ones(sum(L),1) DP_midget_data(L,1)]);
figure; axes; hold on;
plot(DP_midget_data(L,1),DP_midget_data(L,2)/60,'ko');
x = linspace(0,20,10);
plot(x,10.^([ones(length(x),1) x']*b));
ylabel('DF diameter (deg)');
xlabel('Eccentricity (deg)');
title('Dacey and Petersen midget DF (Halve to get 1 SD?)');
set(gca,'Ylim',[0 .2],'Xlim',[0 18])

% % If I assume that the dendritic fields measured by Dacey are 2x a single
% % SD then I get...
% b = regress(log10(DP_midget_data(L,2)/60/2),[ones(sum(L),1) DP_midget_data(L,1)]);
% figure; axes; hold on;
% plot(DP_midget_data(L,1),DP_midget_data(L,2)/60/2,'ko');
% x = linspace(0,20,10);
% plot(x,10.^([ones(length(x),1) x']*b));
% ylabel('DF diameter (1 SD in deg)');
% xlabel('Eccentricity (deg)');
% title('Dacey and Petersen midget');

% Capturing a few points along the fits in Figure 3B
% to estimate how much larger human midget RFs are than 
% monkey midget RFs.
macaque_fit = [3.7834063	1.7960986
5.111495	2.080055
7.5863857	2.7221355
10.003713	3.7911108
12.1781435	4.963961
13.988484	5.9620514
14.9539995	6.5740294];

human_fit = [4.151036	2.135953
5.4205494	2.5678375
6.6294656	3.0492535
7.657095	3.533175
10.013824	4.8603997
11.220212	5.424024
12.608898	6.3594
13.875884	7.184816
15.566041	8.631251];

figure; axes; hold on;
plot(macaque_fit(:,1),log10(macaque_fit(:,2)),'ko');
plot(human_fit(:,1),log10(human_fit(:,2)),'bo');
set(gca,'Yscale','linear'); ylabel('log10 DF diam (min arc)'); xlabel('Ecc. (deg)');
b_macaque = regress(macaque_fit(:,2),[ones(size(macaque_fit,1),1) macaque_fit(:,1)]); % two independent regressions
b_human = regress(human_fit(:,2),[ones(size(human_fit,1),1) human_fit(:,1)]); % two independent regressions
b_macaque(2)/b_human(2) % estimate 1 of scale factor

ecc = [human_fit(:,1); macaque_fit(:,1)];
I_monkey = [zeros(size(human_fit,1),1);ones(size(macaque_fit,1),1)];
X = [ones(size(human_fit,1)+size(macaque_fit,1),1) ecc I_monkey]; % Assuming shift in log DF diam = scaling of DF size
b = regress(log10([human_fit(:,2); macaque_fit(:,2)]),X); % One big regression
plot([2 16]',[1 2 1; 1 16 1]*b,'k-');
plot([2 16]',[1 2 0; 1 16 0]*b,'b-');
10.^b(3)


%%
% Section 6: Derrington and Lennie 1984 (parvocellular)
DL_midget_nasal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984ParvoNasalData'); % x = eccentricity (deg), y = RF field radius (deg)
DL_midget_temporal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984ParvoTemporalData'); % x = eccentricity (deg), y = RF field radius (deg)

figure; axes; hold on;
plot(DL_midget_nasal_data(:,1),DL_midget_nasal_data(:,2),'bo');
plot(DL_midget_temporal_data(:,1),DL_midget_temporal_data(:,2),'r+');
set(gca,'Yscale','log');
b_nasal = regress(log10(DL_midget_nasal_data(:,2)),[ones(length(DL_midget_nasal_data),1) DL_midget_nasal_data(:,1)]);
x = linspace(0,20,100);
plot(x,10.^([ones(length(x),1) x']*b_nasal),'b-');
b_temporal = regress(log10(DL_midget_temporal_data(:,2)),[ones(length(DL_midget_temporal_data),1) DL_midget_temporal_data(:,1)]);
plot(x,10.^([ones(length(x),1) x']*b_temporal),'r-');
xlabel('Eccentricity  (deg)');
ylabel('r_c  (deg)');
title('Derrington and Lennie. parvocellular. Compare to Figure 6A and 6B');

% "*0.61" to convert nasal RFs to equivalent temporal eccentricity
% "*2" to convert radius to diameter (diameter is 2 SDs)
% "/sqrt(2)" to go out only 1 SD instead of 1/e

figure; axes; hold on;
plot(DL_midget_nasal_data(:,1)*.61,DL_midget_nasal_data(:,2)*2/sqrt(2),'bo');
plot(DL_midget_temporal_data(:,1),DL_midget_temporal_data(:,2)*2/sqrt(2),'r+');
set(gca,'Yscale','log');
y = log10([DL_midget_nasal_data(:,2); DL_midget_temporal_data(:,2)].*2./sqrt(2));
X = [ones(length(DL_midget_nasal_data)+length(DL_midget_temporal_data),1) [DL_midget_nasal_data(:,1)*.61; DL_midget_temporal_data(:,1)]];
b_equiv_temporal = regress(y,X);
plot(x,10.^([ones(length(x),1) x']*b_equiv_temporal),'r-');
title('Derrington and Lennie 1984 parvocellular 1 SD');
xlabel('Eccentricity temporal equivalent (deg)');
ylabel('Center diameter (1 SD in deg)');
set(gca,'Yscale','linear');
set(gca,'Xlim',[0 18],'Ylim',[0 .2]);
%%
% Section 7 Derrington and Lennie 1984 (magnocellular)
DL_magno_nasal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984MagnoNasalData'); % x = eccentricity (deg), y = RF field radius (deg)
DL_magno_temporal_data = xlsread('RGC_RFsizes_vs_eccentricity','DL1984MagnoTemporalData'); % x = eccentricity (deg), y = RF field radius (deg)

figure; axes; hold on;
plot(DL_magno_nasal_data(:,1),DL_magno_nasal_data(:,2),'bo');
plot(DL_magno_temporal_data(:,1),DL_magno_temporal_data(:,2),'r+');
set(gca,'Yscale','log');
b_nasal = regress(log10(DL_magno_nasal_data(:,2)),[ones(length(DL_magno_nasal_data),1) DL_magno_nasal_data(:,1)]);
x = linspace(0,20,100);
plot(x,10.^([ones(length(x),1) x']*b_nasal),'b-');
b_temporal = regress(log10(DL_magno_temporal_data(:,2)),[ones(length(DL_magno_temporal_data),1) DL_magno_temporal_data(:,1)]);
plot(x,10.^([ones(length(x),1) x']*b_temporal),'r-');
xlabel('Eccentricity  (deg)');
ylabel('r_c  (deg)');
title('Derrington and Lennie Magnocellular Compare to Figure 6C and 6D');

% "*0.61" to convert nasal RFs to equivalent temporal eccentricity
% "*2" to convert radius to diameter (diameter is 2 SDs)
% "/sqrt(2)" to go out only 1 SD instead of 1/e

figure; axes; hold on;
plot(DL_magno_nasal_data(:,1)*.61,DL_magno_nasal_data(:,2)/sqrt(2),'bo');
plot(DL_magno_temporal_data(:,1),DL_magno_temporal_data(:,2)/sqrt(2),'r+');
set(gca,'Yscale','log');
y = log10([DL_magno_nasal_data(:,2); DL_magno_temporal_data(:,2)]./sqrt(2));
X = [ones(length(DL_magno_nasal_data)+length(DL_magno_temporal_data),1) [DL_magno_nasal_data(:,1)*.61; DL_magno_temporal_data(:,1)]];
b_equiv_temporal = regress(y,X);
plot(x,10.^([ones(length(x),1) x']*b_equiv_temporal),'r-');
title('Derrington and Lennie 1984 magnocellular 1 SD');
xlabel('Eccentricity temporal equivalent (deg)');
ylabel('Center diameter (1 SD in deg)');
set(gca,'Yscale','linear');
set(gca,'Xlim',[0 18],'Ylim',[0 .35]);

% Comparing to Watson model

%%
% Section 8
% Watson 2014. Using parameters from the temporal retina (nasal visual field. 
% Their convention is to represent everything in terms of visual field).
% Midget ganglion cells only (Equation 8)
% remember, this represents both ON and OFF cells toegther. Need to halve
% the density.
a = 0.9851;
r2 = 1.058;
re = 22.14;
dc_0 = 14804.6;
rm = 41.03;
f_0 = 1/1.12;
x = logspace(-1,2,100);
y = 2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)); % From Equation 8
figure; axes; hold on;
loglog(x,y);
set(gca,'Xscale','log','Yscale','log');
title('Watson 2014 midget model. Compare to Figure 9');
xlabel('Eccentricity (deg)');
ylabel('Density (deg^-^2)');
% trying to back out RF size from RF density.
% y is in RF density = RFs/deg^2
% sqrt(1./y(1))

rf_r_deg = logspace(-1,log10(70),100);
a = 0.9729; % Table 1
r2 = 1.084; % Table 1
re = 7.633; % Table 1
dc_0 = 14804.6; % Cone density of fovea
rm = 41.03; % See Equation 7
y = 2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)); % From Equation 8


% Using Equation 9. Distance between adjacent midget RF centers. (assumed to be 2SDs)
rfsize = sqrt(2./(sqrt(3).*y./2)); % Dividing y by 2 to get density of ON (or OFF) cells only
rfsize = rfsize./2; % Dividing whole thing by 2 do get 1 SD

figure; axes; hold on;
plot(rf_r_deg,rfsize);
title('Watson 2014 midget model (1 SD)');
xlabel('Eccentricity temporal equivalent (deg)');
ylabel('RF diameter (1 SD in deg = half inter RF spacing)');
set(gca,'Xlim',[0 20],'Ylim',[0 .2]);

%%
% Section 9: Watanabe and Rodieck 1989 parasol cell DFs
% Having to multiply nasal equivalent eccentricity by 0.61 to get
% everything into *temporal* equivalent eccentricity.

WR_parasol_data = xlsread ('RGC_RFsizes_vs_eccentricity','WatanabeRodieck1989ParasolData'); % x = eccentricity (mm), y = dentritic field diam. (microns)
WR_parasol_data(:,1) = WR_parasol_data(:,1) * 0.61;
figure; axes; hold on;
plot(WR_parasol_data(:,1),WR_parasol_data(:,2),'o');
xlabel('Eccentricity (mm)');
ylabel('Dendritic field diameter (microns)');
title('Watanabe and Rodieck 1989 parasol Compare to Figure 7, bottom, large symbols');
set(gca,'Xlim',[0 6],'Ylim',[0 500],'Ytick',[0:100:500]);

figure; axes; hold on;
plot(WR_parasol_data(:,1)*DEGPERMM,WR_parasol_data(:,2)/1000*DEGPERMM,'o');
b = regress(WR_parasol_data(:,2)*DEGPERMM/1000,[ones(size(WR_parasol_data,1),1) WR_parasol_data(:,1)*DEGPERMM]);
x = linspace(0,18,100);
plot(x,[ones(length(x),1) x']*b);
ylabel('RF diameter (deg)'); % Note this is the full diameter, not 1 SD
xlabel('RF eccentricity (deg)');
title('Watanabe and Rodieck parasol full diameter. Halve to get 1 SD.?');
set(gca,'Xlim',[0 18],'Ylim',[0 0.35])

%%
% Section 9: Watanabe and Rodieck 1989 midget cell DFs

WR_midget_data = xlsread ('RGC_RFsizes_vs_eccentricity','WatanabeRodieck1989MidgetData'); % x = eccentricity (mm), y = dentritic field diam. (microns)
WR_midget_data(:,1) = WR_midget_data(:,1)*0.61;

figure; axes; hold on;
plot(WR_midget_data(:,1),WR_midget_data(:,2),'o');
xlabel('Eccentricity (mm)');
ylabel('Dendritic field diameter (microns)');
title('Watanabe and Rodieck 1989 midget Compare to Figure 7, bottom, small symbols');
set(gca,'Xlim',[0 6],'Ylim',[0 500],'Ytick',[0:100:500]);

figure; axes; hold on;
plot(WR_midget_data(:,1)*DEGPERMM,WR_midget_data(:,2)/1000*DEGPERMM,'o');
b = regress(log10(WR_midget_data(:,2)*DEGPERMM/1000),[ones(size(WR_midget_data,1),1) WR_midget_data(:,1)*DEGPERMM]);
x = linspace(0,18,100);
plot(x,10.^([ones(length(x),1) x']*b));
ylabel('RF diameter (deg)'); % Note this is the full diameter, not 1 SD
xlabel('RF eccentricity (deg)');
title('Watanabe and Rodieck midget full diameter. Halve to get 1 SD.?');
set(gca,'Xlim',[0 18],'Ylim',[0 0.2])

%%
% Section 11
% Cone inverse densities

temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % (Goodchild et al., 1996)
conedensfun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
    (coeffs(3).*(exp(coeffs(4).*x)))+...
    (coeffs(5).*(exp(coeffs(6).*x)));
x = logspace(-1,2,100);
coneDensity = conedensfun(temporalcoeffs,x);
conesPerMM2=coneDensity*1e3;
conesPerDeg = conesPerMM2*MMPERDEG.^2;
figure; axes; hold on;
plot(x,conesPerDeg)
set(gca,'XScale','log','Yscale','log');

rfsize = (pi*sqrt(3))/6.*sqrt(1./(conedensfun(temporalcoeffs,x)*1e3*MMPERDEG.^2));
figure; axes; hold on;
plot(x,rfsize)
set(gca,'XScale','linear','Yscale','linear');
set(gca,'Xlim',[0 18],'Ylim',[0 .1]);

% Wow. The Watson model actually does a very nice job of capturing how cone
% RF diameter changes with eccentricity.