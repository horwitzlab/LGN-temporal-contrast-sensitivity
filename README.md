## Temporal information loss in the macaque visual system


This repository contains data and analysis code from a study of temporal contrast sensitivity. Signal-
noise was compared across four stages of the visual system: photon absorption in cones, cone currents, 
LGN spikes, and behavior. All stimuli used in the LGN recording experiments were at the subjects' 
detection thresholds as described in Gelfand, E. C., & Horwitz, G. D. (2018). Model of parafoveal 
chromatic and  luminance temporal contrast sensitivity of humans and monkeys. Journal of Vision, 18
(12), 1-1; DOI:10.1167/18.12.1.

The data are contained in `LGN_data_stro.mat`. This is a 2x2 Matlab cell array. Two two rows contain 
data from two subjects. The first column contains data from magnocellular neurons, and the second 
column contains data from parvocellular neurons. Within each of the cells in the 2x2 cell array are a 
series of cells, each containing the data from a single neuron.

  * `RUNME.m` is a brief script illustrating calls to `IsoSampGetDPrime` and `IsoSampGetPopulationScaleFactor`.

  * `IsoSampGetDPrime` computes the signal-to-noise of single LGN neurons based on the modulation of the spike rate at the fundamental frequency of the stimulus.

  * `IsoSampGetPopulationScaleFactor` computes a scale factor that approximates how much greater the signal-to-noise of a population of LGN neurons is expected to be than a single one. The single neuron is assumed to be the one recorded and to have a receptive field at the center of the stimulus (consistent with how data were collected). The other cells are assumed to have receptive fields that are arranged on a hexagonal lattice. Signal-to-noise in neurons is assumed to be linear with contrast in the receptive field.
  
  * `IsoSampGetPhotonDPrime` computes the signal-to-noise of an ideal observer of photon catches in the cone photoreceptors.
  
  * `IsoSampConeCurrents` computes the signal-to-noise of an ideal observer of currents across cone photoreceptor cell membranes. Lifted almost verbatim from DTcones.m by Charlie Hass but modified by Greg Horwitz to analyze responses to randomly interleaved temporal frequencies and limit the temporal integration window of the ideal observer.
  
  * `lensDensity_smj.m` Tabulated lens optical density values from Stockman, MacLeod, and Johnson (1993) http://www.cvrl.org/database/text/maclens/lenssmj.htm. Used by IsoSampConeCurrents.m to estimate absorption of light in the lens.
  
  * `MakeFigure5.m` Code that re-creates the main figure of the manuscript, Figure 5, which shows the signal-to-noise ratio at four levels of the visual system: cone photon absorptions, cone currents, LGN populations, and behavioral choices.
  
  * `RGC_RF_sizes.m` Code that compares several estimates of magnocellular (parasol) and parvocellular (midget) RF sizes from the literature. Data are in the Excel spreadsheet `RGC_RFsizes_vs_eccentricity`.
  
 The cone current model requires functions from the Psychophysics Toolbox (http://psychtoolbox.org/).
