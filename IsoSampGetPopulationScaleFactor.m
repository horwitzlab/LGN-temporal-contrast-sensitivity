function [population_scalefactor, ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION,NUMBEROFEYES)
% [population_scalefactor, ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR, RFTRUNCATIONINSD, ONOFFCORRELATION,NUMBEROFEYES)
%
% Find the population scale factor for an ideal observer of LGN spikes.
%
% INPUTS
%
%   stro: structure containing the raw data from a single LGN recording
%   experiment.
%
%   ecc_to_diam_deg: a handle to a function that converts RF eccentricity
%   to diameter (defined as two standard deviations of an assumed Gaussian
%   RF profile).
% 
%   TEMPORONASALSCALEFACTOR: The adjustment to the horizontal RF position
%   to account for the fact that RFs are smaller in the nasal than then
%   temporal retina, RFs sizes returned by ecc_to_diam_deg are in
%   equivalent temporal eccentricity, and we do not know in which retinal 
%   hemifield the RFs of the recorded LGN neurons were located.
%
%   RFTRUNCATIONINSD: How many standard deviations of the Gabor stimulus
%   were displayed in the LGN recording experiments.
%
%   ONOFFCORRELATION: The correlation assumed to exist between ON and OFF
%   LGN mosaics.
%
%   NUMBEROFEYES: 1 for a monocular ideal observer, 2 for binocular. 
%
% OUTPUTS
%
%   population_scalefactor: the factor by which SNR of a population of LGN 
%   neurons is predicting to be greater than the single LGN neuron whose RF
%   is closest to the center of the Gabor stimulus.
%
%   ncells: the number of LGN neurons (in a single mosaic) predicted to have
%   been driven by the Gabor stimulus.

DEBUGGING = 0;

if nargin < 3
    TEMPORONASALSCALEFACTOR = 0.8;
end
if nargin < 4
    RFTRUNCATIONINSD = 4;
end
if nargin < 5
    ONOFFCORRELATION = 0.05;
end
if nargin < 6
    NUMBEROFEYES = 2;
end

bpdf_vec=@(x,y,mu_x,mu_y,sigma)(exp(-((x-mu_x).^2)./(2*sigma^2)-((y-mu_y).^2)/(2*sigma^2))./(2*pi*sigma^2)); % bivariate normpdf
rfx = stro.sum.exptParams.rf_x/10;
rfy = stro.sum.exptParams.rf_y/10;
rf_r_deg = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);
% Cone density is higher in nasal than temporal retina.
% This means that you have to go further out along the nasal retina to get
% to the same cone density as along the temporal retina. e.g. going out 5°
% along the nasal retina gets to the same cone density as going 5*0.61 = 3°
% along the temporal retina. I'm representing RF locations in "equivalent
% *temporal* retinal eccentricity". So I want to divide by 0.61 to get to
% equivalent *nasal* retinal eccentricity (or slightly close to 1 because I
% don't know which side of the retina the RF is on).

RF_diam_deg = ecc_to_diam_deg(rf_r_deg); % 2 SDs of the Gaussian RF
RFdistance = RF_diam_deg; % RF centers are 2 SD apart
RF_STD = RF_diam_deg/2; % 1 standard deviation of Gaussian RF
sigma_gabor = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma'))); % DVA
sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));

x_deg = linspace(-sigma_gabor*2*sigmas_n,sigma_gabor*2*sigmas_n,40);
Rad3Over2 = sqrt(3)/2;
x_centers = [x_deg(1,1):RFdistance:x_deg(1,end)+RFdistance*Rad3Over2];
ncells = length(x_centers)^2;
closest_to_zero = find(abs(x_centers) == min(abs(x_centers)),1);
x_centers = x_centers-x_centers(closest_to_zero);
x_centers_mat = repmat(x_centers,size(x_centers,2),1);
y_centers = x_centers;
x_centers_mat = x_centers_mat*Rad3Over2;
y_centers_mat = repmat(y_centers',1,size(y_centers,2));
y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end) = y_centers_mat(:,rem(find(y_centers == 0),2)+1:2:end)+.5*RFdistance;
mus = zeros(numel(x_centers_mat),1); % Integrated contrast inside truncated Gaussian
interRFdistances = zeros(numel(x_centers_mat),numel(x_centers_mat));

for j = 1:numel(x_centers_mat)
    overlap_point = @(x,y) bpdf_vec(x,y,0,0,RF_STD).*bpdf_vec(x,y,x_centers_mat(j),y_centers_mat(j),sigma_gabor);
    mus(j)=integral2(overlap_point,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD); % RF is truncated at 4 SDs
    
    for k = 1:numel(x_centers_mat) % Tabulating distances between RF centers which we'll need for S2
        interRFdistances(j,k) = sqrt((x_centers_mat(j)-x_centers_mat(k)).^2+(y_centers_mat(j)-y_centers_mat(k)).^2);
    end
end

% filling in S2
disp(['Computing ',num2str(numel(x_centers_mat)),' covariances. Hold please.'])
S2 = nan(numel(x_centers_mat));
for dist = unique(interRFdistances)'
    overlap_point = @(x,y) bpdf_vec(x,y,0,0,RF_STD).*bpdf_vec(x,y,0,dist,RF_STD);
    if (dist >= RFTRUNCATIONINSD*RF_STD)
        S2(interRFdistances==dist)=0;
    else
        S2(interRFdistances==dist)=integral2(overlap_point,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD,RFTRUNCATIONINSD*RF_STD,-RFTRUNCATIONINSD*RF_STD);
    end
end
S2 = S2./max(S2(:)); % Cov depends on doproduct between RF and if two RFs are identical cov = 1 (StatsStuff.m Section 13)

mus = mus./max(mus);
% Normalization ensures that d' of neuron with central RF has the
% observed d'.
%S2 = eye(size(S2)); % DEBUGGING
weights = S2\mus; % these are the ideal weights. No need to normalize to max(weights)
if DEBUGGING
    weights = ones(size(weights));
end
% weights = S2\mus/sqrt(mus'*inv(S2)*mus)
single_mosaic_mean = mus'*weights;
single_mosaic_var=weights'*S2*weights;
population_scalefactor = (2*single_mosaic_mean)/sqrt(2*single_mosaic_var+2*abs(ONOFFCORRELATION)*single_mosaic_var)*sqrt(NUMBEROFEYES);
% Above, in the denominator, I'm adding all four elements of a 2x2
% covariance matrix to get the variance of the sum.
% If ONOFFCORRELATION = 1, including second mosiac doesn't decrease
% the noise at all. If ONOFFCORRELATION = 0, including the second mosaic
% decreases the noise by sqrt(2). In the above equation, I'm *adding* ON and OFF 
% contributions (signal and variances add), taking the mean of the two eyes
% (signal stays the same, but variance decreases by sqrt(2)).


% Brute force sanity check
%invpsf = @(w)sqrt(w'*S2*w)./(mus'*w) % inverse SNR (NSR) to minimize
%[out,fval,exitflag] = fminsearch(invpsf,weights+normrnd(0,.0001,length(weights),1))

% Some debugging code to confirm that my estimates of RF density are the
% same as Watsons 2014 after going through all the numerical manipulations.
% area_of_square_region = (x_deg(1)+x_deg(end)).^2; % Size of square in which RFs and stimulus are placed (underestimate)
%area_of_square_region = (x_centers_mat(end,end)+x_centers_mat(1,end)).^2; % Tight square, overestimate


% % DEBUGGING stuff below
% density = numel(x_centers_mat)./area_of_square_region; % # RFs/deg^2
% %rf_r_deg % eccentricity
% % Watson model
% a = 0.9729; % Table 1
% r2 = 1.084; % Table 1
% re = 7.633; % Table 1
% dc_0 = 14804.6; % Cone density of fovea
% rm = 41.03; % See Equation 7
% y = 2*dc_0.*(1+rf_r_deg./rm).^-1.*(a*(1+(rf_r_deg./r2)).^-2+(1-a)*exp(-rf_r_deg./re)); % From Equation 8
% sprintf('My density estimate: %d, Watson estimate: %d. Ratio: %d',density,y/2, density/(y/2))

end