% Code for computing the signal-to-noise ratio of individial and
% populations of LGN neurons

% Loading the data
load('LGN_data_stro');
stro = data{1}{1}; % Data from a single magnocellular LGN neuron

% IsoSampGetDPrime
% Takes the data from a single neuron (a single stro structure)
% and calculates d' as a function of stimulus condition (color and temporal
% frequency)
[uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro);
mn = nanmean(dprime);
disp(['The SNR for this single neuron is ', num2str(mn)]);

% IsoSampGetPopulationScaleFactor
% Takes the data from a single neuron (a single stro structure) and
% caculates the population scale factor (how much signal-to-noise of a
% population of similar neurons with RFs tile the stimulus is expected 
% to be).

ecc_to_diam_deg = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % magnocellular RF size as a function of eccentricity (temporal retina equivalent)
[population_scalefactor, ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg);

disp(['The population scalefactor is ', num2str(population_scalefactor),' so the populaton SNR is ',num2str(population_scalefactor*mn)])
