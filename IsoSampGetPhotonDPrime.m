function dprime = IsoSampGetPhotonDPrime (flashTimeProfile, frameRate, bkgndlms_Rstar, sigmaInPix, n_cones, uniquestim)
% dprime = IsoSampGetPhotonDPrime(flashTimeProfile, framerate, bkgndlms_Rstar, sigmaInPix, n_cones, uniquestim)
%
% INPUT
%  flashTimeProfile: a vector of contrasts describing the temporal profile
%    of the contrast envelope. By default, the ideal observer integrates over
%    the entire stimulus duration. To limit temporal integration, reduce
%    the length of flashTimeProfile.
%  frameRate: The frame rate of the monitor (1/the time between elements
%    in flashTimeProfile in seconds).
%  bkgndlms_Rstar: number of photons absorbed per L- and M-cone per second in
%    response to the background. A 2-element vector (L and M cones). (Computed by
%    initMonitorCalibration, a subfunction of DTcones.m)
%  sigmaInPix: Standard deviation of the Gabor stimulus in pixels.
%  ncones: nxnx2 matrix of the number of L and M cones per pixel. Depends on
%    eccentricity. (Computed by initConeMosaic, a subfunction of DTcones.m)
%  uniquestim: n x 3 matrix of L-cone contrast, M-cone contrast, temporal frequency.
% 
% OUTPUT
%  dprime: Signal-to-noise ratio of ideal observer with one element per row in uniquestim.
%
% To ensure compatibility between the cone current ideal observer use the
% same values of bkgndlms_Rstar, sigmaInPix, and ncones as used in that
% model. 

nsteps = length(flashTimeProfile);
sperstep = 1./frameRate;
t = linspace(0,nsteps*sperstep,nsteps);

npix = size(n_cones,1);
halfSize = npix/2;
row = -halfSize:halfSize-1;
col = -halfSize:halfSize-1;
row = row-mean(row); 
col = col-mean(col);
[X, Y] = meshgrid(row, col);

dprime = zeros(size(uniquestim,1),1);

for j = 1:size(uniquestim,1)
    sinusoid = sin(2*pi*uniquestim(j,3)*t);
    temporalwtframes = flashTimeProfile.*sinusoid;
    stim2D = exp(-(X.^2 + Y.^2) ./ (2.* sigmaInPix.^2));
    spacetimewt = repmat(stim2D,[1 1 nsteps]).*repmat(permute(temporalwtframes,[1 3 2]),[npix npix 1]);
     
    % Noise
    LMScc = [0 0 0]; % cone contrast at peak
    L = bkgndlms_Rstar(1)*(1+LMScc(1).*spacetimewt).*repmat(squeeze(n_cones(:,:,1)),1,1,nsteps)*sperstep; % L isomerizations each pixel, each frame
    M = bkgndlms_Rstar(2)*(1+LMScc(2).*spacetimewt).*repmat(squeeze(n_cones(:,:,2)),1,1,nsteps)*sperstep; % M isomerizations each pixel, each frame
    Lmn = sum(sum(sum(L.*spacetimewt)));
    Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
    Mmn = sum(sum(sum(M.*spacetimewt)));
    Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2))));
    noisedist = [Lmn Mmn; Lsd Msd];
    
    % Signal
    LMScc = [uniquestim(j,1) uniquestim(j,2) 0]; % cone contrast at peak
    L = bkgndlms_Rstar(1)*(1+LMScc(1).*spacetimewt).*repmat(squeeze(n_cones(:,:,1)),1,1,nsteps)*sperstep; % L isomerizations each pixel, each frame
    M = bkgndlms_Rstar(2)*(1+LMScc(2).*spacetimewt).*repmat(squeeze(n_cones(:,:,2)),1,1,nsteps)*sperstep; % M isomerizations each pixel, each frame
    Lmn = sum(sum(sum(L.*spacetimewt)));
    Lsd = sqrt(sum(sum(sum(L.*spacetimewt.^2)))); % Poisson assumption: var(spacetimewt*X(lambda)) = spacetimewt^2*lambda
    Mmn = sum(sum(sum(M.*spacetimewt)));
    Msd = sqrt(sum(sum(sum(M.*spacetimewt.^2))));
    signaldist = [Lmn Mmn; Lsd Msd];
    
    dprimes = (signaldist(1,:)-noisedist(1,:))./signaldist(2,:);
    dprime(j) = sqrt(dprimes*dprimes');
    
end

end

% Typical values for input parameters are:
% bkgndlms_Rstar: [8898 7378] (L, M)
% sigmaInPix: 4.17 (=0.15 sigma/deg * 27.8 pixels/degree)
% ncones: depends on stimulus size. We assume the L and M cone densities are 
% matched. For a .04 degree Gaussian 10 degrees on 
% the horizontal meridian ncones should be:

% n_cones(:,:,1) % Number of L cones per pixel
%    1.0992    1.1029    1.1068    1.1106
%    1.1005    1.1043    1.1082    1.1120
%    1.1018    1.1057    1.1095    1.1134
%    1.1031    1.1070    1.1109    1.1148

% n_cones(:,:,2) % Number of M cones per pixel
%    1.0992    1.1029    1.1068    1.1106
%    1.1005    1.1043    1.1082    1.1120
%    1.1018    1.1057    1.1095    1.1134
%    1.1031    1.1070    1.1109    1.1148