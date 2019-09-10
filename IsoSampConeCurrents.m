function [gab, cones, mon, idlob, params] = IsoSampConeCurrents(params,saveflag)
%
%  [gab, cones, mon, idlob, params] = IsoSampConeCurrents(params,[saveflag])
%
% 'params' should have the following fields:
%
% params.runType            => 'DTNT', 'DTV1', 'absThresh', 'default'
% params.obsMethod          => ideal observer method (e.g., 'obsMethod_all')
% params.impulseResponse    => 'rieke', or 'deltafxn' 
% params.DTV1_fname     	=> name of a DT expt on V1 (.nex file)
% params.DTNT_fname     	=> name of a DTNT batch file
% params.monCalFile     	=> the name of a monitor calibration file (e.g., DTcones.mat)
% params.unitTest       	=> true or false (runs the testing routines)
% params.eqMosaic       	=> true or false (equates the cone mosaics)
% params.saveDir            => path to the directory where data will be saved
% params.notes          	=> a string. notes for the data file that describes the experimental run
% params.parallelOperations => true/false. so that DTcones knows how to store temporary files used during parfor loops
%
% Created by CAH using cone physiology data (and fits) from Juan Angueyra



    %create/initialize the structures that will define the simulation
    [gab, params] = initGabor(params);
    eyes = initEyeParams(params);
    [mon, params] = initMonitorCalibration(params, gab, eyes);
    cones = initConeMosaic(gab, mon, eyes);
    [cones, gab] = makeConeLinearFilter(params, cones, mon, gab);
    cones = defineConeGainScaleFactor(cones, mon);
    cones = makeConePowerSpectrum(cones, gab, params);
    mon = getBkgndLinearResp(cones, mon);
    idlob = initIdealObserver(params.obsMethod, gab, cones, mon); %this is where you set which ideal observer will be used

    % for absThres experiments, eliminate most of the cones, and set
    % the bkgnd Rstar to zero.
    if strcmpi(params.runType, 'absthresh')
        [cones, mon] = setupAbsThreshParams(cones, mon, params); 
    end 
    
    % for comparisons to V1 neurophysiology, one option is to modify the
    % cone mosaic to fit the RF size of the V1 neuron. Restrict the cone
    % mosaic here:
    if isfield(params, 'aperatureMosaic') && params.aperatureMosaic
        cones = aperatureConeMosaic(params, cones, mon);
    end
    
   
    %stop here and perform some testing routines if desired
    if params.unitTest
        callAllTestingSubfunctions(gab, mon, cones, params, idlob);
        return %exit DTcones
    elseif params.eqMosaic
        [~,cones] = feval(idlob.method, 'equate_mosaic', idlob, gab, cones);
    end

    %loop through each color direction specified. Make a gabor stimulus and
    %determine the output of the cone mosaic. Do this in subfunctions so
    %that the temporary variables created do not clog up the workspace (and
    %eat RAM).
   
    
    nColors = size(gab.colorDirs, 1);
    for clr = 1:nColors
        fprintf('  Color %d of %d\n', clr, nColors);
        nContrasts = length(gab.contrasts{clr});

        for cnt = 1:nContrasts

            fprintf('   Contrast %d of %d\n', cnt, nContrasts);
            
            % create a gabor movie
            lms_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
            gab.driftRate =  gab.driftRates{clr}(cnt);
            gab.movie = makeGaborMovie(lms_Rstar, gab, mon, 'movie', cones, params, 0); %stimulus movie in 4D [X, Y, Time, Cones], where Cones=3. Should be in units of R*/sec due to the stimulus

            % calculate the linear response of the cones
            [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar);
            

            % sythesize multiple trials by adding noise to the linear
            % response and then performing the ideal observer analysis.
            % This is the olddb (Monte-Carlo) style. If gab.nTrials is set to
            % zero, then no Monte-Carlo trials will be calculated and only
            % the analytic solution will be derived...
            for trl = 1:gab.nTrials
                if ~rem(trl, 100); fprintf('Trial %d \n', trl); end
                
                cones = shapedConeNoise(cones, gab);
                [idlob, cones] = feval(idlob.method, 'compute', idlob, gab, cones, mon, clr, cnt, trl, params); %computes the ideal observer response
            end

            % Now compute the variance of the mosaic's population response
            % using the analytic method
            [idlob, cones] = feval(idlob.method, 'analytic_solution', idlob, gab, cones, mon, clr, cnt, trl, params);
        end
    end

    %save the results of the experiment in a directory. This directory
    %should contain the results, all relavant simulation structures, AND an
    %exact copy of the code that was run to produce the results. Time stamp
    %the directory name. Re-name the .m file, or add a different extension
    %so that it is identifiable as the actual code ran.
    if saveflag
        saveDataToDisk(gab, cones, mon, idlob, params)
    end

end



%
%   SUPPORTING FUNCTIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gab, params] = initGabor(params)
    
    % common to both methods:
    nContrasts = 22;
    gab.nTrials = 0;
    
    switch lower(params.runType)
        
        case 'default'
            gab.sd = 3;           % in 1/10 dva
            gab.nSd = 4;          % number of SDs in gabor
            gab.theta = 0;        % in radians. Zero is horizontal drifting up
            gab.gamma = 1;        % SDx == SDy
            gab.sf = 1;           % in cpd
            gab.driftRate = 3;    % in cyc/sec
            gab.length = 0.666;   % in sec
            gab.rf_x = -50;        % in 1/10 dva
            gab.rf_y = -35;         % in 1/10 dva
            gab.colorDirs = [1, 1, 1;...
                1, -1, 0;...
                0, 0, 1];
            cntrsts = [0, logspace(log10(0.0001), log10(0.15), nContrasts)];
            gab.contrasts = repmat({cntrsts}, size(gab.colorDirs,1), 1);
            
            
        case 'dtv1'
            
            % unpack the .nex file
            DT = dtobj(params.DTV1_fname);
            DT = stripOutGratingTrials(DT);
            if DT.sum.exptParams.expt_meth ~= 1;
                error(' The specified file is not from a MoCS experiment')
            end
            
            %assign the relavant parameters
            gab.sd = unique(DT.trial(:, DT.idx.gaborSigma));
            gab.nSd = DT.sum.exptParams.flash_size;
            gab.theta = unique(DT.trial(:, DT.idx.gaborTheta));
            gab.gamma = unique(DT.trial(:, DT.idx.gaborGamma));
            gab.sf = 1 ./ (unique(DT.trial(:, DT.idx.gaborLambda)) ./ DT.sum.exptParams.pixperdeg);
            gab.driftRate = unique(DT.trial(:, DT.idx.driftRate));
            gab.length = DT.sum.exptParams.flash_length./1000;
            gab.rf_x = DT.sum.exptParams.rf_x;
            gab.rf_y = DT.sum.exptParams.rf_y;
            
            %determine the color directions used
            gab.colorDirs = reshape(DT.sum.exptParams.RF_colors, [], 3)';
            blankIdx = sum(abs(gab.colorDirs), 2) == 0;
            gab.colorDirs(blankIdx,:) = [];
            norms = sqrt(sum(gab.colorDirs.^2,2));
            gab.colorDirs = bsxfun(@rdivide, gab.colorDirs, norms); %now unit vecs
            gab.colorDirs(gab.colorDirs(:,1)<0,:) = -gab.colorDirs(gab.colorDirs(:,1)<0,:); %flip the signs for L-cone weights that were negative
            l_siso = softEq([0,0,1], abs(gab.colorDirs), [], 'rows');
            if sum(gab.colorDirs(l_siso,:)) < 0;
                gab.colorDirs(l_siso,:) = gab.colorDirs(l_siso,:) .* -1;
            end
            
            
            %now determine the rgb's used
            x = 0:255; %the normal range of the gamma look up table
            xx = linspace(0, 255, 2^16); %the desired quantization of the gammaTable
            g1 = reshape(DT.sum.exptParams.gamma_table, 256, 3);
            gammaTable = [spline(x, g1(:,1), xx)', spline(x, g1(:,2), xx)', spline(x, g1(:,3), xx)'];
            M = reshape(DT.sum.exptParams.m_mtx,3,3);
            bkgndrgb = [DT.sum.exptParams.bkgnd_r, DT.sum.exptParams.bkgnd_g, DT.sum.exptParams.bkgnd_b];
            bkgndlms = M * bkgndrgb(:);
            for clr  = unique(DT.trial(:, DT.idx.colorDir))'
                for cnt = unique(DT.trial(:, DT.idx.cntrstLev))'
                    if cnt == 1
                        tList = DT.trial(:, DT.idx.cntrstLev) == cnt;
                    else
                        tList = (DT.trial(:, DT.idx.colorDir) == clr) & (DT.trial(:, DT.idx.cntrstLev) == cnt);
                    end
                    
                    %pull out the RGB for this clr/contrast
                    RGB = unique(DT.trial(tList, [DT.idx.flashR, DT.idx.flashG, DT.idx.flashB]), 'rows');
                    if size(RGB, 1) ~= 1; error('more than one RGB per color direction'); end
                    rgb = [gammaTable(RGB(1), 1), gammaTable(RGB(2), 2), gammaTable(RGB(3), 3)];
                    lms = M * rgb(:);
                    V1_CCs{clr}(:,cnt) = norm((lms-bkgndlms)./bkgndlms);
                end
                
                % instead of sampling the same contrasts as the DT expt. Resample a
                % range from very low contrast to the max contrast used in the expt
                gab.contrasts{clr} = [0, logspace(log10(0.00001), log10(V1_CCs{clr}(end)), nContrasts)];
            end
            
            % package the stro file (here named 'DT') so that other
            % functions can have access to it
            params.V1stro = DT;
            
        case 'dtnt'
            
            % load the batch dtnt data
            load(params.DTNT_fname);
            in = dtnt;
            
            % specify all the gabor parameters.
            gab.rf_x = in.rf_x;
            gab.rf_y = in.rf_y;
            gab.sd = in.sigma;
            gab.nSd = in.nSD;
            gab.theta = in.theta;
            gab.gamma = in.gamma;
            gab.length = in.length;
            gab.driftRate = in.speed;
            
            gab.colorDirs = in.colorDirs;
            gab.sf = in.sfs;
            gab.contrasts = {};
            for a = 1:numel(in.alphas)
                gab.contrasts{a} = [0, logspace(log10(0.0003), log10(in.alphas(a)), nContrasts)];
            end
            
        case 'absthresh'
            
            % Most of these parameters don't matter b/c later I'll zero out
            % the entire stiulus except for a single pixel. The only thing
            % that doesn't get changed is the stimulus length. Also, the
            % contrast no longer signifies CC, and now signifies an integer
            % number of R* (to be delivered in a single time step)
            gab.sd = 1;           % in 1/10 dva
            gab.nSd = .2;         % number of SDs in gabor
            gab.theta = 0;        % in radians. Zero is horizontal drifting up
            gab.gamma = 1;        % SDx == SDy
            gab.sf = 1;           % in cpd
            gab.driftRate = 3;    % in cyc/sec
            gab.length = 0.400;   % in sec. should be same length as the filter
            gab.rf_x = 25;        % in 1/10 dva
            gab.rf_y = 25;        % in 1/10 dva
            gab.colorDirs = [1, 1, 1];
            gab.contrasts = {0:1:275}; % integer number of R*, NOT cone contrast!!!
      
        case 'isosamp'
            if ~isfield(params,'stro')
                error('cannot find stro structure in params stucture');
            end
            sigma_idx = strcmp(params.stro.sum.trialFields(1,:),'sigma');
            sigmas_n_idx = strcmp(params.stro.sum.trialFields(1,:),'sigmas_n');
            sf_idx = strcmp(params.stro.sum.trialFields(1,:),'sf');
            tf_idx = strcmp(params.stro.sum.trialFields(1,:),'tf');
            flashtime_idx = strcmp(params.stro.sum.trialFields(1,:),'flash_time');
            theta_idx = strcmp(params.stro.sum.trialFields(1,:),'theta');
            gamma_idx = strcmp(params.stro.sum.trialFields(1,:),'gamma');
            stim_l_idx = strcmp(params.stro.sum.trialFields(1,:),'stim_l');
            stim_m_idx = strcmp(params.stro.sum.trialFields(1,:),'stim_m');
            
            if ~isfield(params, 'gab') % Using the std of the stimulus
                gab.sd = unique(params.stro.trial(:,sigma_idx))*10;  % Charlie's code expects gabor sd to be in 1/10 dva
                gab.nSd = unique(params.stro.trial(:,sigmas_n_idx));
            else 
                gab.sd = params.gab.sd*10; % Interpreting params.gab.sd as the sd of the Gaussian field of cones we're integrating over (could be RF of a single neuron)
                gab.nSd = 2;
            end
            gab.theta = unique(params.stro.trial(:,theta_idx));
            gab.gamma = unique(params.stro.trial(:,gamma_idx));
            gab.sf = unique(params.stro.trial(:,sf_idx));
            TF = params.stro.trial(:,tf_idx);
            Lcc = params.stro.trial(:,stim_l_idx);
            Mcc = params.stro.trial(:,stim_m_idx);
            uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
            
            cdir = uniquetol(atan2(uniquestim(:,2),uniquestim(:,1)));
            cdir(cdir == 0) = [];
            gab.contrasts = {};
            gab.driftRates = {}; % Hack: driftRates are all the TFs and driftRate (no "s")
            gab.driftRate = 0; % will be updated in main function before the call to makeGaborMovie.
          
            for i = 1:size(cdir,1)
                L = softEq(atan2(uniquestim(:,2),uniquestim(:,1)), cdir(i));
                % Below three lines are unnecessary. Noise mean is zero.
                %if i ==1 & any(all(uniquestim == 0,2)) == 1 % assigning zero contrast condition to cdir 1
                %    L = L | all(uniquestim == 0,2);
                %end
                gab.contrasts{length(gab.contrasts)+1} = sqrt(sum(uniquestim(L,[1 2]).^2,2));
                gab.driftRates{length(gab.driftRates)+1} = uniquestim(L,3);
            end
           
            gab.colorDirs = [cos(cdir) sin(cdir) zeros(size(cdir,1),1)]; % L M S
            gab.length = unique(params.stro.trial(:,flashtime_idx))/1000; % in sec
            
            % 8/9/19 passing in params.temporalIntegrationLims
            if ~isfield(params,'temporalIntegrationLims') | isempty(params.temporalIntegrationLims)
                params.temporalIntegrationLims = [0 gab.length];
            end
            if length(params.temporalIntegrationLims) ~= 2
                error('params.temporalIntegrationLims must be a two-element vector');
            end
            gab.rf_x = params.stro.sum.exptParams.rf_x;
            gab.rf_y = params.stro.sum.exptParams.rf_y;
    end
end

function movie = makeGaborMovie(LMS_Rstar, gab, mon, TYPE, cones, params, initPhase)

    % Define the constants
    nPix = round((gab.nSd*2) .* (gab.sd./10) .* mon.pixperdeg);

    halfSize = round(nPix./2);
    gaborPhase = initPhase;
    flashDeltaPhase = gab.driftRate * 2 * pi * (1./mon.frameRate); %the amount to advance each frame
    sigmaInPix = (gab.sd./10) .* mon.pixperdeg;
    lambdaInPix = 1./ (gab.sf ./ mon.pixperdeg);
    
    % Parameters for ramping on/off stimulus
    % GDLH 7/23/18: flashNumFrames is the duration of the full stimulus
    % (ramp on, plateau, ramp off). Adjusting flashTimeProf to be the
    % temporal profile of the stimulus we want to use in the ideal
    % observer simulation (which could be longer or shorter than the
    % full stimulus).

    % added 7/23/18 for backward compatability GDLH
    if ~isfield(gab,'flashTimeProf')
        flashNumFrames = ceil(mon.frameRate * gab.length);
        rampLength = ceil(flashNumFrames / 4); % in frames
        ramp = linspace(0, 1, rampLength);  % ramp is 1/4th of the total duration on either side
        plateau = ones(1,flashNumFrames - (rampLength .* 2));
        gab.flashTimeProf = [ramp, plateau, fliplr(ramp)];
    end
    flashNumFrames = length(gab.flashTimeProf);
    % end of 7.23.18 edits
    
    % meshgrid for the gabor
    row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window; 
    col = -halfSize:halfSize-1;
    row = row-mean(row); % 8/5/19. Centering Gabor GDLH 
    col = col-mean(col);
    [X, Y] = meshgrid(row, col);
    xprime = X .* cos(-gab.theta) + Y .* sin(-gab.theta);
    yprime = -X .* sin(-gab.theta) + Y .* cos(-gab.theta);

    % Preallocate space for the movie (4D: [X, Y, Time, Guns])
    movie = nan(size(X,1), size(X,1), flashNumFrames, 3);
    
    % This function can be used to make a stimulus movie, or a
    % spatiotemporal weighting function for ideal observer analysis, or a
    % simple stimulus for the absolute threshold expts.
    if strcmpi(TYPE, 'movie') && ~strcmpi(params.runType, 'absthresh') % the normal case, just making a gabor
            Rstar_increment = LMS_Rstar - mon.bkgndlms_Rstar;
            bkgnd_Rstar = mon.bkgndlms_Rstar;
    elseif strcmpi(TYPE, 'idlob template')
            bkgnd_Rstar = [0 0 0];
            Rstar_increment = LMS_Rstar; % should be 1, -1, or 0.
    elseif strcmpi(TYPE, 'movie') && strcmpi(params.runType, 'absthresh') % for absThresh expts
            movie(:) = 0;
            movie(1,1,1,1) = LMS_Rstar(1); % put an integer number of R* in the first pix, first time step, L cone map.
            return % nothing else to do...
    end
    
    
    % Make the movie frame by frame.
    for t = 1:flashNumFrames

        % Multiply the gabor by the increment from background and by the temporal weighting fxn.
        gabor = exp(-(xprime.^2 + gab.gamma.^2 .* (yprime).^2) ./ (2.* sigmaInPix.^2)) .* cos(2 .* pi .* (yprime) ./ lambdaInPix + gaborPhase);

        %R*/sec for each pix/frame/cone
        movie(:,:,t,1) = (gabor .* gab.flashTimeProf(t) .* Rstar_increment(1)) + bkgnd_Rstar(1);
        movie(:,:,t,2) = (gabor .* gab.flashTimeProf(t) .* Rstar_increment(2)) + bkgnd_Rstar(2);
        movie(:,:,t,3) = (gabor .* gab.flashTimeProf(t) .* Rstar_increment(3)) + bkgnd_Rstar(3);

        % Update the phase
        gaborPhase = gaborPhase + flashDeltaPhase;
    end
    
    % adjust the ideal observer template for the visual latency. Add a
    % bunch of zeros to the front and then take off an equal number on the
    % back.
    if strcmpi(TYPE, 'idlob template') && ...
            any(strcmpi(params.obsMethod, {'obsMethod_noClrEqSpace', 'obsMethod_absThresh'}));
        disp('latency corrected')
        latencyInFrames = round(cones.latency.*mon.frameRate);
        latency_movie = zeros([size(movie,1), size(movie,2), latencyInFrames, 3]);
        movie = cat(3, latency_movie, movie);
        movie(:,:,flashNumFrames+1:end,:) = [];
    end
end

function lms_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones)
    % this function will calculate the LMS_Rstar values for each
    % color/contrast combination. If the color dirs are specified as
    % strings, than assume the user wants gun iso colors
    
    if ischar(gab.colorDirs(clr,:))
        
        % gun iso specified as 'rgun', 'ggun', 'bgun' in the 'gab.colorDirs
        % field.
        switch gab.colorDirs(clr,:)
            case 'rgun'
                gunisodir = [1 0 0];
            case 'ggun'
                gunisodir = [0 1 0];
            case 'bgun'
                gunisodir = [0 0 1];
        end
        
        % change the gunisodir into LMS_Rstar.
        GC = gab.contrasts{clr}(cnt); % Gun Contrast
        gunDirContrast = gunisodir(:) .* GC;
        rgbTrial = mon.bkgndrgb(:) .* (gunDirContrast + 1);
        lms_Rstar = mon.rgb2Rstar * rgbTrial; % convert from gun intensity space to cone R* space
        
    elseif strcmpi(params.runType, 'absthresh')
        % in this case, the "contrast" array specifies an integer number of
        % photons...
        if cones.samplingRate ~= 1000 || mon.frameRate ~= 100;
            mon.frameRate
            cones.samplingRate
            error('Sampling rate or refresh rate not set appropriately')
        end
        pulseWidth = 0.010; %only works if refresh and cone sampling rates are set correctly
        RStarPerSec = gab.contrasts{clr}(cnt) ./ pulseWidth; %an integer number of photons delivered in one time step (expressed in R*/sec)
        lms_Rstar = [RStarPerSec; 0; 0]; % only deliver the stimulus to the Lcone mosaic
    
    else % 'default', 'dtnt', or 'dtv1' type runs...
        
        % the color dirs are specified in LMS CC units. This corresonds to
        % a vector in LMS CC space, where the space is defined by what ever
        % cone fundamentals are specified by DTcals.mat. To emulate the
        % ACTUAL stimuli used in a physiology/behvioral experiment, I'm
        % goint to take the LMS values specified by the user (in
        % gab.colorDirs) and determine what rgb directions this would be
        % for the monitor used in the monkey experiments. Then I'll convert
        % into LMS_R* units for the cone model. This way, the cone model is
        % tested on the identical stimuli (in rgb units) as the monkeys
        
        % determine the LMS color dir asked for and force it to be a unit
        % vector
        colorDir_generic = gab.colorDirs(clr,:);
        colorDir_generic = colorDir_generic ./ norm(colorDir_generic); %make sure it's a unit vector
        
        % determine the LMS cone contrast color and convert this into a
        % vector that specifies the lms values of the stimulus
        CC = gab.contrasts{clr}(cnt);
        colorDir_CC = colorDir_generic .* CC;
        bkgndlms_deviceSpecific_smj = mon.Mmtx * mon.bkgndrgb(:);
        gaborlms_deviceSpecific_smj = bkgndlms_deviceSpecific_smj .* (1+colorDir_CC(:));
        
        % determine what rgb direction this would correspond to on the
        % moniotor used in the monkey experiments
        gaborrgb_deviceSpecific = mon.Mmtx \ gaborlms_deviceSpecific_smj;
        
        % convert the device specific rgbs into lms units for the cone
        % model
        lms_Rstar = mon.rgb2Rstar * gaborrgb_deviceSpecific(:);
        
        % a simple error check on the cone sampling rate
        if rem(cones.samplingRate,mon.frameRate) > 1e-2 % cones.samplingRate should be (nearly) an integer multiple of the mon refresh rate
            warning('might have the incorrect cone sampling rate')
        end
        
        % make sure that none of the lms_Rstar values are negative.
        diff_from_bkgnd = lms_Rstar(:)-mon.bkgndlms_Rstar(:);
        assert(~any(diff_from_bkgnd > mon.bkgndlms_Rstar(:)), sprintf('ERROR: R* values went negative for color <%d> cnt <%d>', clr, cnt))
    end

end



function eyes = initEyeParams(params)
    switch lower(params.eyeType)
        case 'monkey'
            eyes.diam = 19;                 % in mm; From Qiao-Grider 2007: (19mm)
            eyes.pupilArea = 12.566;        % in mm^2... assuming a radius of 2mm which is pretty close to Freya's measurements (12.566)
            eyes.postFocalLength = 12.75;   % in mm. From Qiao-Grider 2007 (12.75)
            eyes.fovealMacpig = 0.35;        % Wooten 2005 says 0.3, but Snoderly suggests it's higher... (default = 0.35)
            
        case 'human'
            error('No human values set')
    end
    if isfield(params,'eyeNumber')
        eyes.number = params.eyeNumber;
    else
        eyes.number = 1;
    end
end

function [mon, params] = initMonitorCalibration(params, gab, eyes)

    % load in some information about the monitor.
    switch lower(params.runType)
        case {'default', 'dtnt', 'absthresh','isosamp'}
            %load in calibration data from a representative file
            if isstruct(params.monCalFile)
                cal = params.monCalFile;
                fprintf(' Calibration structure in use: %s \n', cal.fname)
            else
                fprintf(' Calibration file in use: %s \n', params.monCalFile)
                load(params.monCalFile)
            end
            % pacakge some other relavant info
            monSpect = reshape(cal.monSpect, [],3);
            mon.bkgndrgb = cal.bkgndrgb;
            mon.frameRate = cal.frameRate;
            mon.monSpectWavelengths = 380:5:780;
            mon.pixperdeg = cal.pixperdeg;
            mon.Mmtx = reshape(cal.Mmtx,3,3);
            
            % in the case where you load in equalBkgnd cal file, make sure
            % that the stimulus RF position is what we expect. Due to
            % macpig, equalBkgnd must be specified uniquely for each
            % spatial location
            if isfield(cal, 'rfX')
                assert(cal.rfX==gab.rf_x, '  ->->->  ERROR: for equal bkgnd, RF needs to be at (-5 -3.5) DVA')
                assert(cal.rfY==gab.rf_y, '  ->->->  ERROR: for equal bkgnd, RF needs to be at (-5 -3.5) DVA')
            end
            
        case 'dtv1'
            DT = params.V1stro;
            monSpect = reshape(DT.sum.exptParams.mon_spect, [], 3);
            mon.bkgndrgb = [DT.sum.exptParams.bkgnd_r, DT.sum.exptParams.bkgnd_g, DT.sum.exptParams.bkgnd_b];
            mon.frameRate = DT.sum.exptParams.frame_rate;
            mon.monSpectWavelengths = 380:5:780;
            mon.pixperdeg = DT.sum.exptParams.pixperdeg;
            mon.Mmtx = reshape(DT.sum.exptParams.m_mtx,3,3);
            params = rmfield(params, 'V1stro'); % no longer need the V1 stro struct
    end
    
    % the monitor spectra are in units of spectral radiance (W/(str*M^2)).
    % Convert these values to retinal irradiance in W/um^2.
    monitorIrradiance = retinalIrradiance(monSpect, eyes.pupilArea, eyes.diam);
    
    % do some pre-retinal filtering on the monitor (irradiance spectra).
    % Filter through the lens, and through some macular pigment
    monitorIrradiance = preRetinalFilter_lens(monitorIrradiance);
    monitorIrradiance = preRetinalFilter_macpig(monitorIrradiance, gab, mon, eyes);
    
    % "monitorIrradiance" is in units of W/um^2/sec. Now convert to
    % photons/um^2/sec. The energy of a photon is equal to:
    %
    %   Energy of photon = hc/lambda
    %
    % where h = plank's constant, c = speed of light, lambda = wavelength
    % in meters.
    %
    % To convert to photons/s/cone we just divide the energy due to the
    % monitor by the energy in a single photon:
    %
    %   photons = monitorIrradiance ./ (hc/lambda);
    hc = 6.626069e-34 .* 3e8; %Planck's constant * Speed of Light
    lambdas = repmat(mon.monSpectWavelengths', 1, 3) ./ 1e9; %converting from nm to meters
    photonSpectra = monitorIrradiance ./ (hc./lambdas);
    
    
    %generate the cone action spectra using the wavelengths at which the
    %monitor was calibrated
    actionSpectra(1,:) = coneActionSpectra('LCone',mon.monSpectWavelengths); % in normalized units???
    actionSpectra(2,:) = coneActionSpectra('MCone',mon.monSpectWavelengths);
    actionSpectra(3,:) = coneActionSpectra('SCone',mon.monSpectWavelengths);
    
    % convert to absorPTance, and incorporate the cone collecting area.
    % Adding the cone collecting area after the normalization means that
    % the fundamentals peak at the value of the cone collecting area.
    opticaldensity = 0.3; % From SMJ (1993) and other papers. Bowmaker 1978
    fundamentals = 1 - 10.^(-actionSpectra .* opticaldensity);
    fundamentals = bsxfun(@rdivide, fundamentals, max(fundamentals,[],2));
    coneCollectingArea = 0.6; %in um^2 (JA & CAH).. about 0.6 Schneeweis (1999) and Schnapf (1990)
    
    fundamentals = fundamentals .* coneCollectingArea;
    
    % make a matrix that converts from gun (intensity) space to cone (R*/sec)
    % space:
    mon.rgb2Rstar = fundamentals * photonSpectra;
    
    % determine the R* due to the background intensity of the monitor
    mon.bkgndlms_Rstar = mon.rgb2Rstar * mon.bkgndrgb(:);
    
    % generate a rod-absorPTance spectum and determine the rod R* due to the background
    rodActionSpectra = coneActionSpectra('rod',mon.monSpectWavelengths);
    rod_OD =  0.35; % Rodieck estimate approx 0.47, Bowmaker et al (1978) estimates 0.475, Baylor (1984) 0.3 to 0.4;
    rod_fund = 1 - 10.^(-rodActionSpectra .* rod_OD);
    rod_fund = bsxfun(@rdivide, rod_fund, max(rod_fund,[],2));
    rodCollectingArea = 1; %in um^2 (JA & CAH).. about 1 Schneeweis (1999) and Schnapf (1990)
    rod_fund = rod_fund .* rodCollectingArea;
    
    rgb2RodRstar = rod_fund * photonSpectra;
    mon.bkgndrod_Rstar = rgb2RodRstar * mon.bkgndrgb(:);
    
end

function irradiance = retinalIrradiance(radiance,pupilAreaMM,eyeSizeMM)
    %
    % Perform the geometric calculations necessary to convert a measurement of source
    % radiance to corresponding retinal irradiance. 
    %
    %   Input radiance should be in units of power/m^2-sr-wlinterval.
    %   Input pupilAreaMM should be in units of mm^2.
    %   Input eyeSizeMM should be the length of the eye in mm.
    %   Output irradiance is in units of power/um^2-sec-wlinterval.
    %
    % Jan 2012 (CAH) Stealing code from David Brainard. See ptb's function
    %                "RadianceToRetIrradiance" for more details.
    
    % Convert power/sr-M^2-wlinterval to power/sr-mm^2-wlinterval
    radianceMM = radiance.*1e-6;
    
    % There are 4pi steradians per Surface Area of a sphere.
    % Since SA = 4*pi*r^2, 1 stradian = r^2. If the eye is circumscribed by
    % a hemisphere, than the number of steradians that the pupil covers is
    % simply the pupilArea/(eyeDiam^2);
    steradian = eyeSizeMM^2; %in mm^2
    nSteradians = pupilAreaMM/steradian;
    irradianceMM = radianceMM .* nSteradians;

    % Convert units to um^2 from mm^2 base.
    irradiance = irradianceMM.*1e-6;
 end

function out = preRetinalFilter_lens(in)
    
    % define how the lens density should be manipulated
    METHOD = 'horwitz'; % 'Horwitz' for Gregs estimates based on monkey psychophysics, or 'SMJ' for the typical human measurements
    
    %load in the lens density measurements
    lensDensity_smj
    
    switch lower(METHOD)
        case 'horwitz'
            %make the lens density at 400nm equal to one. Arbitrary, but
            %results in good fits to the monkey psychophysics
            idx = find(S_lens_smj == 400);
            newDensity = den_lens_smj./den_lens_smj(idx);
            lenstransmittance = 1./(10.^(newDensity));
        case 'smj'
            % divide the lens density measurements by 1.16 to convert to
            % "open pupil", and multiply by 1.28 as per SMJ 1993
            scaleFactor = 1.28 ./ 1.16;
            newDensity = den_lens_smj .* scaleFactor;
            lenstransmittance = 1./(10.^(newDensity));
    end
    
    % filter the monSpect by multiplying by the lens transmitance 
    out = bsxfun(@times, in, lenstransmittance);
    
end

function newIrradiance = preRetinalFilter_macpig(oldIrradiance, gab, mon, eyes)
    
    % calculate the gabor's eccentricity
    gaborEccentricity = sqrt((gab.rf_x./10)^2 + (gab.rf_y./10)^2);
    
    % determine the pigment density according to published data in humans
    % (Wooten 2005) with modifications for lower density in monkeys.
    macpigDensity = eyes.fovealMacpig .* exp(-gaborEccentricity./1.03); % from Wooten 2005, Snodderly (for macpig in monkeys)
    
    % load in the macular pigment absorbance spectrum. Scale the spectrum
    % to reflect the macular pigment density at the gabor eccentricity.
    % Next, convert to tranmittance.
    load 'den_mac_ws' % measured at the same wavelengths as the monitor...
    macpig_wf = S_mac_ws(1):S_mac_ws(2):(S_mac_ws(1)+S_mac_ws(2)*S_mac_ws(3)-1);
    if ~all(macpig_wf == mon.monSpectWavelengths); error('Mac pig wavelengths mismatch'); end;
    idx_460 = macpig_wf == 460;
    macpigAbsorbance = den_mac_ws ./ den_mac_ws(idx_460); % normalized to height = 1 at 460 nm
    macpigAbsorbance = macpigAbsorbance .* macpigDensity; % set's the absorbance at 460 nm according to the gabor eccentricity
    macpigTransmittance = 1./(10.^macpigAbsorbance); % transmittance is the reciprocal of absorbance
    
    % modify the monitor irradiance to reflect the macular pigment
    newIrradiance = bsxfun(@times, oldIrradiance, macpigTransmittance);
        
end

function Spectrum=coneActionSpectra(Type,Wavelength)
    %
    % Spectrum=coneActionSpectra(Type,Wavelength)
    %
    % Type options: 'LCone','MCone','SCone','Rod'
    %
    % Generates Normalized Photoreceptor Spectral Sensitivity Curve relying on
    % fits derived from Baylor's et al. suction recordings (1985, 1987) on
    % macaque photoreceptors and fitted according to the Govardovskii nomogram
    % (2000). Fitting procedure can be checked in FittingGovadorvskii.m
    %
    % Dec_2010 (Angueyra) Created the function.
    % Jan 2012 (CAH) Incorporated Juan's code into this script under the name coneActionSpectra

    switch lower(Type)
        case 'lcone'
            %values derived from fitting procedure
            LambdaMax=565.2836;
            BetaBandScaling=0.1460;
            Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);
        case 'mcone'
            %values derived from fitting procedure
            LambdaMax=534.3201;
            BetaBandScaling=0.1797;
            Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);        
        case 'scone'
            %values derived from fitting procedure
            LambdaMax=430.4053;
            BetaBandScaling=0.5210;
            Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);  
        case 'rod'
            %values derived from fitting procedure
            LambdaMax=493.3022;
            BetaBandScaling=0.1233;
            Spectrum=GovadorvskiiTemplateA1([LambdaMax,BetaBandScaling],Wavelength);
    end
end

function fit = GovadorvskiiTemplateA1(beta, wavelength)

    % Govardovskii template for A1-based pigments with alpha nad beta band
    % Based on Govadorvskii et al., 2000
    %
    % Dec 2010 (Angueyra) Created the code
    % Jan 2012 (CAH) Pasted the code into this script without modification

    A=69.7;
    a=0.88;
    B=28.0;
    b=0.922;
    C=-14.9;
    c=1.104;
    D=0.674;

    alphafit=1./(exp(A.*(a-(beta(1)./wavelength)))+exp(B.*(b-(beta(1)./wavelength)))+exp(C.*(c-(beta(1)./wavelength)))+D);


    lambdaMaxBeta=189+0.315*beta(1);
    b=-40.5+0.195*beta(1); %bandwidth
    betafit=beta(2)*exp(-(((wavelength-lambdaMaxBeta)./b).^2));

    fit=alphafit+betafit;
end

function cones = initConeMosaic(gab, mon, eyes)
    

    gaborWidth_pix = round((gab.nSd*2) .* (gab.sd./10) .* mon.pixperdeg);
    
    if rem(gaborWidth_pix,2)
        gaborWidth_pix = gaborWidth_pix+1; % make this an even number for consistency with the Gabor movie
    end
    gaborWidth_dva = gaborWidth_pix./mon.pixperdeg;
    pixWidth_dva= gaborWidth_dva ./ gaborWidth_pix;
    
    % Assume that a principle ray from the bottom of the pixel and a ray
    % from the top of the pixel each travel through the nodal point of the
    % lens. This means that the DVA of the pixel will be equal to the angle
    % subtended by the retinal image (measured at the nodal point).
    W_ret = tand(pixWidth_dva) .* eyes.postFocalLength;
    pixRetinalArea_mm2 = W_ret^2;
    
    % now find the eccentricity of each of the pixel centers, and solve for
    % the cone density at each of these locations. then allocate the
    % correct number of cones to these locations based off the retinal area
    % of each pixel.
    gaborHalfWidth_dva = gaborWidth_dva ./ 2;
    center_x = gab.rf_x/10; %in dva
    center_y = gab.rf_y/10; %in dva
    pix_x = linspace(center_x-gaborHalfWidth_dva, center_x+gaborHalfWidth_dva, gaborWidth_pix*2+1);%finding half pixels
    pix_x = pix_x(2:2:end); %Pixel centers in dva;
    pix_y = linspace(center_y-gaborHalfWidth_dva, center_y+gaborHalfWidth_dva, gaborWidth_pix*2+1);%finding half pixels
    pix_y = pix_y(2:2:end); %Pixel centers in dva;
    
    %make a meshgrid to find each pixel's eccentricity
    [X, Y] = meshgrid(pix_x, pix_y);
    %eccInDVA = sqrt(X.^2 + Y.^2);
    
    %Solving for cone density per pixel. Also taken from juan's ConeGrid function.
    %coeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % for macaque (Goodchild et al., 1996)
    %coneDensity= (coeffs(1).*(exp(coeffs(2).*eccInDVA)))+...
    %    (coeffs(3).*(exp(coeffs(4).*eccInDVA)))+...
    %    (coeffs(5).*(exp(coeffs(6).*eccInDVA)));
    %conesPerMM2=coneDensity*1e3; % in cones/mm^2 Left out in Goodchild et al., 1996? Compare to their Fig 5 & 7
    
    % Assuming that the cone density along the horzontal meridian is
    % the average of the nasal and temporal densities and that the cone
    % density along the vertical merdian is the temporal density. GDLH
    % 9/26/17
    temporalcoeffs=[150.9676 -1.2220 35.9979 -0.1567 9.9936 -0.0258]; % (Goodchild et al., 1996)
    nasalcoeffs=[176.6624 -7.9473 94.4908 -0.3518 18.6761 -0.0236]; % (Packer et al. 1989) See GrantBrainStorming.m section 30
    conedensfun = @(coeffs,x)(coeffs(1).*(exp(coeffs(2).*x)))+...
        (coeffs(3).*(exp(coeffs(4).*x)))+...
        (coeffs(5).*(exp(coeffs(6).*x)));
    [theta,r] = cart2pol(X,Y);
    coneDensity = (cos(theta).^2.*(conedensfun(temporalcoeffs,r)+conedensfun(nasalcoeffs,r))/2) + (sin(theta).^2.*conedensfun(temporalcoeffs,r));
    conesPerMM2=coneDensity*1e3;
    % number of cones per "pixel" on the retina
    conesPerPix_all = conesPerMM2 .* pixRetinalArea_mm2;

    % Solving for S-cone density per MM2. Taken from Monasterio 1985 (only
    % good b/w 2 and 60 dva.
    if any(r(:)<2) || any(r(:)>60)
        warning('Scone density not accurate outside of 2 to 60 dva!')
    end
    mmPerDeg = W_ret .* mon.pixperdeg;
    a = 121.86; % Solving for S-cone density (temporal side).
    b = -0.206;
    c = 89.95;
    d = -0.045;
    SconesPerDeg2 = (a.*exp(b.*r)) + (c.*exp(d.*r));
    SconesPerMM2 = SconesPerDeg2 ./ (mmPerDeg^2);
    conesPerPix_S = SconesPerMM2 .* pixRetinalArea_mm2;
    
    %assemble the outputs
    cones.num_S = conesPerPix_S;
    cones.num_L = (conesPerPix_all - conesPerPix_S) ./ 2;
    cones.num_M = (conesPerPix_all - conesPerPix_S) ./ 2;
    
    if eyes.number == 2
        % give the model 2 eyes
        warning('model has 2 eyes')
        cones.num_L = cones.num_L .* 2;
        cones.num_M = cones.num_M .* 2;
        cones.num_S = cones.num_S .* 2;
    else
        warning('model has 1 eye')
    end
end

function cones = aperatureConeMosaic(params, cones, mon)
    % used only for cones/V1 comparisions. This function limits the size of
    % the cone mosaic to be roughly the size of a V1 RF. This function will
    % not change the dimensionality of the stimulus (which is not
    % computationally effecient) but will allow the rest of the code to run
    % un-modified. I'll just put zeros in the cone mosaic where there
    % should be no cones.
    
    % right now, the aperature module is only available for the filtered
    % weight function ideal observer. I'll need to update the other ideal
    % observer modules following the methodology used in
    % 'obsMethod_filteredWtFxn'.
    assert(strcmpi(params.obsMethod, 'obsMethod_filteredWtFxn'), 'ERROR: aperature setting is currently only available for the ''filtered weight function'' observer');
    
    % load in the GT data file
    GT = gtobj(params.GTV1_fname);
    spikenum = strcmp(GT.sum.rasterCells(1,:),getSpikenum(GT, 'first'));
    cones.aptDiam = getAperatureTuning(GT, spikenum);
    
    % make the aperature. The aperature size is in units of DVA, so work in
    % DVA...
    halfSize = size(cones.num_L, 1)/2; % in pix
    row = (-halfSize:halfSize-1) ./ mon.pixperdeg; % now in DVA
    col = (-halfSize:halfSize-1) ./ mon.pixperdeg;
    [X, Y] = meshgrid(row, col);
    aptur = sqrt((X.^2+Y.^2)); % the distance of each pixel from the center of the stim (in DVA)
    idx = aptur>(cones.aptDiam/2); % extend the aperature one radius out.
    aptur(idx) = 0; % set the flanks to zero.
    aptur(~idx) = 1; % set everything else to 1.
    
    % clear out the cone mosaic outside of the V1 RF.
    cones.num_L = cones.num_L .* aptur;
    cones.num_M = cones.num_M .* aptur;
    cones.num_S = cones.num_S .* aptur;

end

function [cones, gab] = makeConeLinearFilter(params, cones, mon, gab)
    %
    %Stealing code from juan to make the linear filter for the cones. I'll
    %export the filter in the fourier domain b/c this is the only version
    %that gets used in coneVolution.
    
    %figure out the cone "sampling rate"
    cones.samplingRate = params.coneSampRate; %desired rate in Hz. Good candidates: [525 600 675 750 825 900 975]. These all give rise (nearly) to an integer number of cone samples per monitor refresh
    cones.samplesPerFrame = round(cones.samplingRate./mon.frameRate);
    cones.samplingRate = cones.samplesPerFrame.*mon.frameRate; %adjusted due to rounding
    
    % determine how many samples the cone impulse response should contain.
    % Also determine how many samples to add to the beginning of the Gabor
    % to allow the simulation to attain steady state.
    gab.framesInNewGabor = ceil(mon.frameRate * gab.length) .* cones.samplesPerFrame; %up sample due to the cone sampling rate
    filterTimeLength = 0.400; %seconds
    cones.filterSamples = ceil(filterTimeLength .* cones.samplingRate);
    cones.nFrontPad = cones.filterSamples+5; %add to the front of the gabor to allow the system to come to steady state
    cones.convLength = (cones.nFrontPad + gab.framesInNewGabor) + cones.filterSamples - 1; % conv gives length N+M-1
    cones.nBackPad = cones.convLength - (cones.nFrontPad + gab.framesInNewGabor); % add to the back of the GABOR so that the dimensionality of the gabor stimulus is 1xConvLength
    %make the filter, and estimate the visual latency. Juan says the units
    %are time (in sec) vs pA/R*/cone
    TimeAxis=[0:(cones.convLength-1)] ./ cones.samplingRate;
    ScFact = 0.6745;% coef(1);
    TauR = 0.0216; %Rising Phase Time Constant
    TauD = 0.0299; %Damping Time Constant
    TauP = 0.5311; %Period
    Phi = 34.1814; %Phase
    Filter = ScFact .* (((TimeAxis./TauR).^3)./(1+((TimeAxis./TauR).^3))) .* exp(-((TimeAxis./TauD))).*cos(((2.*pi.*TimeAxis)./TauP)+(2*pi*Phi/360));

    % estimate latency
    [~, idx] = max(Filter);
    cones.latency = TimeAxis(idx);
    
    % for the expts where I need to use a delta fxn for the IRF, just put a
    % 'one' where the visual latency typically is. The height of the IRF
    % will get scaled later (as per usual).
    if strcmpi(params.impulseResponse, 'deltafxn')                
            Filter = zeros(1, cones.convLength);
            latencyInFrames = round(cones.latency .* cones.samplingRate);
            Filter(latencyInFrames) = 1; % this gets adjusted by a gain scale factor later (to the height that a normal IRF would be scaled given the current adaptation state)
    end
    
    % many of the samples in the filter are just there as place holders to
    % make the dimensionality consistent with the gabor (to facilitate FFT
    % convolution and mulitplication in the fourier domain). These values
    % should be close to zero, but set them explicitly to zero. Doing so
    % makes the FFT-conv identical to the real conv (see testing routines).
    Filter(cones.filterSamples+1:end) = 0;
    
    % convert to fourier domain
    cones.filter_fft = fft(Filter); %not normalizing here, but will compensate in coneVolution_FFT (below)
end

function cones = defineConeGainScaleFactor(cones, mon)
    %
    % The IRF used by the model was calculated at ~4500 R*. Make sure to
    % adjust the gain of the model's IRF so that the model gain is
    % consistent with their adaptation state
    
    % adjust the gain by assuming cone gain is a Weber-Fechner relationship
    % with the half-desensitizing value indicated in Juan's paper.
    Io = 2250;                 % half-desensitizing value (in R*/sec, from Juan's paper) [old value was 4500]
    Ib = mon.bkgndlms_Rstar;   % from the monitor initalization routine
    gain_dark = 0.32;          % from Juan's paper (approximate, and in units of pA/R*)[old value was 0.16]
    gainRatio = 1 ./ (1+(Ib./Io));
    gainAtBkgnd = gainRatio .* gain_dark;
    gainOfIRF = max(ifft(cones.filter_fft));
    cones.gainScaleFactor = gainAtBkgnd ./ gainOfIRF;    

end

function cones = makeConePowerSpectrum(cones, gab, params)
    % Using Juan's code to generate cone noise based on the power spectrum
    % calculated from recordings of cone noise fit with a sum of 2
    % Lorentzians up to 1kHz.
    %
    % Angueyra  1/2010  Created the function
    % CAH       2/2012  Integrated the fxn into DTcones and changed how the
    %                   frequency axis is defined for odd and even N.
    %           9/2012  Switching to two-sided PS
    %

    % constants:
    lorentzCoeffs = [0.2, 30, 2.0, 0.05, 180, 2.5]; % from Juan's experiments
    
    samplingRate = cones.samplingRate; %in samp/sec

    % Calculate the frequncy axis
    N = gab.framesInNewGabor;
    if rem(N,2)
        k = -((N-1)./2):((N-1)./2); % pos freqs when N is Odd. For two sided freqs: k = -((N-1)./2):((N-1)./2)
    else
        k = -(N./2):((N./2)-1); % 0:(N./2) => abs(neg Freqs) when N is even. For two sided freqs: k = -(N./2):((N./2)-1)
    end
    freqAxis = (k./N).*samplingRate;

    % generate the positive portion of the PS
    l_posFreq = freqAxis>=0;
    tmpFreqs = freqAxis(l_posFreq);
    PS_posFreq = abs(lorentzCoeffs(1)) ./ (1 + (tmpFreqs ./ abs(lorentzCoeffs(2))).^2).^lorentzCoeffs(3);
    PS_posFreq = PS_posFreq + abs(lorentzCoeffs(4)) ./ (1 + (tmpFreqs ./ abs(lorentzCoeffs(5))).^lorentzCoeffs(6));

    % generate the negative portion of the PS by finding the PS at the
    % abs(negative freqs). This is b/c juan's model is only defined on the
    % positve interval
    l_negFreq = freqAxis<0;
    tmpFreqs = abs(freqAxis(l_negFreq));
    PS_negFreq = abs(lorentzCoeffs(1)) ./ (1 + (tmpFreqs ./ abs(lorentzCoeffs(2))).^2).^lorentzCoeffs(3);
    PS_negFreq = PS_negFreq + abs(lorentzCoeffs(4)) ./ (1 + (tmpFreqs ./ abs(lorentzCoeffs(5))).^lorentzCoeffs(6));
    
    
    % Put the neg and pos parts of the PS together. Divide by 2 so that
    % the units are for a two-sided PS. 
    cones.modelNoise_ps = [PS_posFreq, PS_negFreq];
    cones.modelNoise_ps = cones.modelNoise_ps ./ 2;

    % Multiply by delta Frequency so
    % that the units are in pow/bin... not pow/Hz.
    dF = unique(diff(freqAxis));
    dF = dF(1);
    cones.modelNoise_ps = cones.modelNoise_ps .* dF;
       
    % DON'T SET THE DC POWER TO ZERO. Even though there isn't any noise at
    % the DC, becuase the PS is sampled discretely, setting the power at DC
    % to zero makes the power effectively zero for all frequencies in the
    % first bin. If the first bin is from 0 to 2.5 Hz, than these non-zero
    % frequencies will have less noise than would be expected...
    if cones.modelNoise_ps(1) == 0; error('Noise at DC set to zero'); end
    
    
    % When Juan makes the PS I think normalizes-out the effect of discrete
    % sampling. Bring it back so that the values of the PS are similar to
    % those used natively in matlab. In particular, the model's PS should
    % behave such that the sum(PS)/N^2 should have a value equal to the
    % variance in the time domain....
    cones.modelNoise_ps = cones.modelNoise_ps ./ (1./(N.*(N-1))); % norm fact for conversion from PS to sample var
    
    % make the PS flat (if needed), but maintain the Variance (according to parseval's
    % theorem for a discrete process. These lines of code are rather
    % verbose, but I like that at least I know why and how things are
    % getting normalized.
    if params.flatPowerSpect
        oldVar = sum(cones.modelNoise_ps) .* (1./(N*(N-1)));
        cones.modelNoise_ps = ones(size(cones.modelNoise_ps));
        newVar = N .* (1./(N*(N-1)));
        cones.modelNoise_ps = ones(size(cones.modelNoise_ps)) .* (oldVar./newVar);
        sum(cones.modelNoise_ps) .* (1./(N*(N-1)))
        warning('using a flat noise PS')
    end
    
    % turn the power spectra into ordinary magnitudes:
    cones.modelNoise_fft = sqrt(cones.modelNoise_ps);
    cones.freqAxis_fft = [freqAxis(freqAxis>=0) freqAxis(freqAxis<0)];
   
end

function [convMovie, movie] = coneVolution_FFT(movie, gab, cones, bkgndRstar, DEBUG)
    
    % front pad, back pad, and upsample the gabor movie.
    [convMovie, movie] = upsampleAndBookend(movie, cones, bkgndRstar);
    
    if exist('DEBUG', 'var') && strcmpi(DEBUG, 'impulse_response')
        % make an impulse in the middle of the stimlus. (adding this here
        % b/c i need access to a stimulus in 'cone sampling time'
        convMovie(:,:,:,:) = 0;
        idx = round( cones.nFrontPad + round(gab.framesInNewGabor./2) );
        convMovie(:,:,idx,:) = cones.samplingRate;
    end
    
    %convert the gabor movie to the fourier domain. the gabor movie is now
    %contained in convMovie to allow for in-place operations.
    convMovie = fft(convMovie, [], 3); %run the FFT in the "time" dimension
    
    % multiply the gabor and the coneFilter in the fourier domain.
    convMovie = bsxfun(@times, permute(cones.filter_fft, [1,3,2]), convMovie);
    
    % convert back to the time domain and compensate for the sampling rate.
    % The impuseresponseFFT and the stimulusFFT have not been adjusted for
    % the sampling rate, so do it here:
    %
    % ifft(FFT(x)/rate .* FFT(y)/rate)*rate = ifft(FFT(x) .* FFT(y)) ./ rate
    %
    convMovie = ifft(convMovie, [], 3) ./ cones.samplingRate;
    
    % hack off the bookends
    tRange = (cones.nFrontPad+1) : (size(convMovie,3) - cones.nBackPad);
    convMovie = convMovie(:,:, tRange, :);
    
    % adjust for the gain of the cones
    convMovie = bsxfun(@times, convMovie, permute(cones.gainScaleFactor, [4,3,2,1]));
    
    % simple error checking
    if gab.framesInNewGabor ~= size(convMovie, 3); error(' #### Size Mismatch in ConeVolution'); end
    
end

function [upsamp, gaborMovie] = upsampleAndBookend(gaborMovie, cones, bkgndRstar)

    % start by allocating space to an array that will become the new gabor
    % movie. This array should have the proper number of "time" dimensions
    % to allow the system to come to steady state, to acomodate the cone
    % sampling rate, and to avoid circular convolution via the FFT method.
    dims = size(gaborMovie);
    upsamp = nan(dims(1), dims(2), cones.convLength, dims(4)); %this will become the new gabor movie
    upsamp(:, :, 1:cones.nFrontPad, 1) = bkgndRstar(1);
    upsamp(:, :, 1:cones.nFrontPad, 2) = bkgndRstar(2);
    upsamp(:, :, 1:cones.nFrontPad, 3) = bkgndRstar(3);
    
    
    % fill in the gabor movie and upsample to the cone sampling rate
    idx = cones.nFrontPad+1;
    for a = 1:dims(3); % dims(3) = nFrames in old gabor
        frames = idx:(idx+cones.samplesPerFrame-1);
        upsamp(:,:,frames,1) = repmat(gaborMovie(:,:,a,1), [1, 1, cones.samplesPerFrame]);
        upsamp(:,:,frames,2) = repmat(gaborMovie(:,:,a,2), [1, 1, cones.samplesPerFrame]);
        upsamp(:,:,frames,3) = repmat(gaborMovie(:,:,a,3), [1, 1, cones.samplesPerFrame]);
        
        % update the idx
        idx = idx+cones.samplesPerFrame;
    end
    
    
    % pad each timeseries on the back to avoid circular convolution. The
    % padding must be zeros!!!
    for i = 1:3
        upsamp(:, :, idx:cones.convLength, i) = 0;
    end
    
    % assign the newmovie to "upsamp" (happen when this function
    % gets called) so that the FFT convolution can be done "in place". This
    % aviods making two giant arrays...
    gaborMovie = []; %kill the variable to save RAM
    
end

function mon = getBkgndLinearResp(cones, mon)
    
    % generate a time series that is composed of the bkgndlms_Rstar
    bkgnd = repmat(mon.bkgndlms_Rstar, 1, numel(cones.filter_fft));
    bkgnd(:, end-cones.nBackPad:end) = 0; %the padding for circular convolution needs to be zeros...
    bkgnd_fft = fft(bkgnd,[],2);
    linresp = bsxfun(@times, cones.filter_fft, bkgnd_fft) ./ cones.samplingRate;
    linresp = ifft(linresp, [], 2);
    linresp = linresp(:, (cones.nFrontPad+1) : (cones.convLength-cones.nBackPad));
    linresp = mean(linresp, 2); % reducing to 3 numbers and taking out the numerical noise
    
    % comensate for the gain of the cones
    mon.bkgndlms_pA = linresp .* cones.gainScaleFactor;
        
end

function cones = shapedConeNoise(cones, gab, DEBUGMODE)
    
    % FFT of gaussian noise.
    N = gab.framesInNewGabor;
    if N ~= size(cones.linresp, 3); error('size mismatch between noise length and length of linear response'); end
    normfact_ff = 1./(N.*(N-1)); % norm fact to get "sample variance" from parseval's theorem.
    noiseSigma = sqrt(sum(cones.modelNoise_ps).*normfact_ff);
    WN_fft = fft(normrnd(0,noiseSigma,size(cones.linresp)),[],3);
    
    % Shaping White Noise by the PowerSpectrum of Cone Dark Noise
    tmp = permute(cones.modelNoise_fft, [1,3,2]); %so that samples go in the "time" dimension
    cones.noise = bsxfun(@times, tmp, WN_fft);
    
    % Convert back to time domain.
    cones.noise = ifft(cones.noise, [], 3);

    % Normalize Noise to match original variance. Here are some thoughts
    % about the scale factor: 
    % Parseval's theorm for finite length time sequences, as quoted form
    % numerical recipes, "The integral of the one-sided PSD (per unit time)
    % over positve frequency is equal to the mean squared amplitude of the
    % time signal".
    scaleFactor = sqrt( (sum(cones.modelNoise_ps).*normfact_ff) ./ (var(cones.noise, 0, 3)) );
    cones.noise = bsxfun(@times, cones.noise, scaleFactor);
    
    % run some testing routines if need be. Typically DEBUGMODE is only
    % defined if this function is called by a testing subfunction
    if (nargin==3) && strcmpi(DEBUGMODE, 'unit_test')
        runTestSubFxn()
    end
 
    function runTestSubFxn
        
        %
        % COMPARE THE PREDICTED AND AVERAGE NOISE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fft of the noise alone.
        pow = abs(fft(cones.noise, [], 3)); % actually magnitudes
        pow = (pow.^2); % now power
        pow = permute(pow, [3,1,2,4]);
        pow = reshape(pow, N, []); % freqs goes down the columns
        pow_avg = mean(pow, 2)';
        
        
        %plot the average power spectra of the "shaped white noise"
        figure, hold on,
        plot(cones.freqAxis_fft(cones.freqAxis_fft>=0), pow_avg(cones.freqAxis_fft>=0), '-k', 'linewidth', 7)
        plot(cones.freqAxis_fft(cones.freqAxis_fft>=0), cones.modelNoise_ps(cones.freqAxis_fft>=0), 'r-', 'linewidth', 3)
        set(gca, 'yscale', 'log', 'xscale', 'log')
        ylim([850, 85000])
        xlabel('Frequency')
        ylabel('Power per Hz')
        legend('Average noise', 'Model noise')
        drawnow

        
        %
        % MAKE SURE THAT THE AVERAGE POWER SPECT OF WN IS ROUGHLY FLAT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp = permute(WN_fft, [3,1,2,4]);
        tmp = reshape(tmp, N, []); %time goes down the columns
        tmp = abs(fft(tmp,[],1));
        tmp = tmp.^2;
        tmp = mean(tmp,2);
        figure
        set(gcf, 'name', 'Average Powerspectra for White Noise');
        plot(cones.freqAxis_fft(cones.freqAxis_fft>=0), tmp(cones.freqAxis_fft>=0), '-k', 'linewidth', 2)
        ylim([0, max(tmp).*1.4])
        drawnow
        
        %
        % EXAMPLE AUTOCORRELATION FUNCTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nRows = 5;
        pixels = [unidrnd(size(cones.noise,1), nRows^2,1), unidrnd(size(cones.noise,2), nRows^2,1), unidrnd(size(cones.noise,4), nRows^2,1)];
        
        
        figure
        set(gcf, 'position', [148 5 1184 824])
        set(gcf, 'name', 'Auto-correlation of shaped and unshape noise');
        for a = 1:size(pixels,1)
            tmp = permute(cones.noise(pixels(a,1), pixels(a,2), :, pixels(a,3)), [1,3,2,4]);
            [AC_shaped, lags_shaped] = xcorr(tmp, N./2, 'coeff');
            tmp = ifft(permute(WN_fft(pixels(a,1), pixels(a,2), :, pixels(a,3)), [1,3,2,4]));
            [AC_WN, lags_WN] = xcorr(tmp, N./2, 'coeff');
            
            subplot(nRows, nRows, a), hold on
            plot(lags_WN, AC_WN, 'r');
            plot(lags_shaped, AC_shaped, 'k', 'linewidth', 2);
            title(sprintf('[%d, %d, :, %d]', pixels(a,1), pixels(a,2), pixels(a,3)));
            axis tight
        end
        drawnow
        
        %
        % PRINT OUT SOME STATISTICS REGARDING THE SHAPED NOISE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %ttests for shaped noise
        tmp = permute(cones.noise, [3,1,2,4]);
        tmp = reshape(tmp, N,[]);
        [~,p_shaped] = ttest(tmp);
        
        %ttest for WN
        tmp = ifft(WN_fft,[],3);
        tmp = permute(tmp, [3,1,2,4]);
        tmp = reshape(tmp, N, []);
        [~,p_WN] = ttest(tmp);
        
        %descriptive stats
        noiseVar = var(cones.noise, [], 3);
        noiseVar = mean(noiseVar(:));
        noiseMean = mean(cones.noise,3);
        noiseMean = mean(noiseMean(:));
        
        %print the results to the command window
        fprintf('  Some statistics regarding the shaped noise: \n')
        fprintf('  1) Parseval''s theorem assumes zero mean: \n')
        fprintf('     Mean noise is %g pA \n', noiseMean);
        fprintf('     T-tests on shaped noise: p(type I error) = %.3f\n', sum(p_shaped<0.05)./numel(p_shaped));
        fprintf('     T-tests on white noise: p(type I error) = %.3f\n', sum(p_WN<0.05)./numel(p_WN));
        fprintf('  2) Variance of noise is %g pA^2 \n', noiseVar);
        
    end
end



function idlob = initIdealObserver(observerMethod, gab, cones, mon)
    
    fxnHandle = str2func(observerMethod);
    idlob.method = fxnHandle;
    idlob = feval(fxnHandle, 'initialize', idlob, gab, cones, mon); %sets up the 'resp' field of the structure
    
end

function [idlob, cones] = obsMethod_noClrEqSpace(command, idlob, gab, cones, mon, clr, cnt, trl, params)
    %
    % This ideal observer only has access to the spatiotemporal weighting
    % function. Also, the space/time wt function is a gabor stimulus that
    % peaks/troughs at �1. All the cones in the wt fxn are in phase with
    % eachother. The result is a 3 vector for each trial. I'll have to do
    % linear discriminant inorder to determine the neurometric threshold.
    
    switch lower(command)    
        case 'initialize'
            nColors = size(gab.colorDirs,1);
            nContrasts = size(gab.contrasts{1},2);
            nTrials = gab.nTrials;
            
            % dimensionality will be (nColors x nContrasts x nTrials) for
            % the monte carlo method, and (nColors x nContrasts) for the
            % analytic solution. Each entry will be a CELL ARRAY of three
            % numbers, one for each cone type. I'm doing this so that the
            % arrays have the same dimensionality as obsMethod_all and so I
            % can use the same offline analysis scripts.
            [idlob.resp(1:nColors, 1:nContrasts, 1:nTrials)] = deal({nan(1,3)});
            [idlob.analyticMean(1:nColors, 1:nContrasts), idlob.analyticVar(1:nColors, 1:nContrasts)] = deal({nan(1,3)});
            
        case 'compute'
            % subtract off the mean from the linear response
            zeroMeanLinResp = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            
            % generate the pooled response and add the (pooled) noise.
            sigPlusNoise = zeroMeanLinResp + pooledNoise_average(cones.noise, cones.num_L, cones.num_M, cones.num_S);
            clear zeroMeanLinResp %kill to save ram
            
            % weight the sigPlusNoise by a template based on the gabor
            % stimulus. The template will peak at +/- 1 instead
            % of the peak Rstar. The template cone mosaics are all IN PHASE
            % despite the fact that the stimulus may not be.
            switch params.enableScones
                case true
                    wt = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
                case false
                    wt = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, 0);
            end
            [wt, ~] = upsampleAndBookend(wt, cones, [0 0 0]);
            tRange = (cones.nFrontPad+1) : (size(wt,3) - cones.nBackPad);
            wt = wt(:,:, tRange, :); %hack off the bookends.
            
            dotWithLinResp = sum(sigPlusNoise .* wt, 3); %sum across the 'time' dimension
            idealObsResp = sum(sum(dotWithLinResp,1),2); %sum across the spatial dimensions.
            idlob.resp{clr, cnt, trl} = squeeze(idealObsResp);
            
        case 'analytic_solution'
            % compute the mean and variance independently. This will mean
            % incorporating the wt fxns seperately for the mean and variance.
            
            % (1) make a spatiotemporal wt fxn. The template will peak at
            % +/- 1 instead of the peak Rstar.
            switch params.enableScones
                case true
                    wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
                case false
                    wt_spaceTime = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, 0);
            end
            [wt_spaceTime, ~] = upsampleAndBookend(wt_spaceTime, cones, [0 0 0]);
            tRange = (cones.nFrontPad+1) : (size(wt_spaceTime,3) - cones.nBackPad);
            wt_spaceTime = wt_spaceTime(:,:, tRange, :); %hack off the bookends.

            
            
            % (2) calculate the mean response of the ideal observer
            idlOb_mean = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            idlOb_mean = idlOb_mean .* wt_spaceTime;
            idlOb_mean = sum(sum(sum(idlOb_mean,1),2),3); %sum across space. The sum over time is actually part of the dot product with the wt. fxn...
            idlob.analyticMean{clr, cnt} = squeeze(idlOb_mean);
            idlob.wt_spaceTime = wt_spaceTime;

            % (3) calculate the variance of the ideal observer. This
            % should be the same for each contrast (within color
            % directions) so just calculate it on the first trial and then
            % carry it forward for the other contrasts
           % if cnt == 1;
                
                wt_spaceTime_FFT = fft(wt_spaceTime,[],3);
                noise_fft = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_fft, [1,3,2])); % weighting by the space/time/clr fxn in the fft domain
                
                % normalizing by N due to discrete fourier transform (see
                % Plancherel's theorem)
                N = size(wt_spaceTime, 3);
                noise_fft = noise_fft./N;
                
                % finish the 'dot product' in the fourier domain by summing
                % across frequencies
                noise_fft = sum(abs(noise_fft).^2, 3); % variance of a single pixel across time
                noise_fft = squeeze(noise_fft); % remove the (singelton) time dimension
                
                % pool across cones (within a pixel) and then across pixels
                numCones = cat(3, cones.num_L, cones.num_M, cones.num_S);
                noise_fft = noise_fft ./ numCones; % the noise is independent from cone to cone (within pixels) so the average var = var(x)/numCones
                noise_fft = sum(sum(noise_fft,1),2); % ideal observer sums signals from adjecent pix, so variances add.
                
                idlob.analyticVar{clr, cnt} = squeeze(noise_fft);
                
           % else
           %     idlob.analyticVar{clr, cnt} = idlob.analyticVar{clr, 1};
           % end
            
        case 'equate_mosaic'
            
            % For cone isolating stimuli, the ideal observer should produce
            % the identical responses if the cone mosaics are ths same. Use
            % this functionality to run some debugging code.
            cones.num_M = cones.num_L;
            cones.num_S = cones.num_L;
    end
end

function [idlob, cones] = obsMethod_filteredWtFxn(command, idlob, gab, cones, mon, clr, cnt, trl, params)
    %
    % This ideal observer only has access to the spatiotemporal weighting
    % function. Also, the space/time wt function is a gabor stimulus that
    % peaks/troughs at �1. All the cones in the wt fxn are in phase with
    % eachother. The result is a 3 vector for each trial. I'll have to do
    % linear discriminant inorder to determine the neurometric threshold.
    
    switch lower(command)    
        case 'initialize'
            nColors = size(gab.colorDirs,1);
            nContrasts = size(gab.contrasts{1},2);
            nTrials = gab.nTrials;
            
            % dimensionality will be (nColors x nContrasts x nTrials) for
            % the monte carlo method, and (nColors x nContrasts) for the
            % analytic solution. Each entry will be a CELL ARRAY of three
            % numbers, one for each cone type. I'm doing this so that the
            % arrays have the same dimensionality as obsMethod_all and so I
            % can use the same offline analysis scripts.
            [idlob.resp(1:nColors, 1:nContrasts, 1:nTrials)] = deal({nan(1,3)});
            [idlob.analyticMean(1:nColors, 1:nContrasts), idlob.analyticVar(1:nColors, 1:nContrasts)] = deal({nan(1,3)});
            
        case 'compute'
            
            % make sure the user doesn't trust the results when the cone
            % mosaic is aperatured and the user wants monte carlo trials
            % too.
            assert(~params.aperatureMosaic, 'ERROR: Monte carlo trials are not supported when the cone aperature setting is activated')
            
            % subtract off the mean from the linear response
            zeroMeanLinResp = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            
            % generate the pooled response and add the (pooled) noise.
            sigPlusNoise = zeroMeanLinResp + pooledNoise_average(cones.noise, cones.num_L, cones.num_M, cones.num_S);
            clear zeroMeanLinResp %kill to save ram
            
            % weight the sigPlusNoise by a template based on the gabor
            % stimulus. The template will peak at +/- 1 instead
            % of the peak Rstar. The template cone mosaics are all IN PHASE
            % despite the fact that the stimulus may not be.
            switch params.enableScones
                case true
                    wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
                case false
                    wt_spaceTime = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, 0);
            end
            
            % filter the wt fxn through the cone IRF
            [wt_spaceTime, ~] = coneVolution_FFT(wt_spaceTime, gab, cones, [0 0 0]);
            
            % the function that convolves the wtFxn also adjust the output
            % according to the cone gain. This is a hack to undo this
            % normalization so that the wt fxn peaks/troughs at �1.
            scaleFact = max(max(max(abs(wt_spaceTime),[],1),[],2),[],3);
            wt_spaceTime = bsxfun(@rdivide, wt_spaceTime, scaleFact);
            
            
            dotWithLinResp = sum(sigPlusNoise .* wt_spaceTime, 3); %sum across the 'time' dimension
            idealObsResp = sum(sum(dotWithLinResp,1),2); %sum across the spatial dimensions.
            idlob.resp{clr, cnt, trl} = squeeze(idealObsResp);
            
        case 'analytic_solution'
            % compute the mean and variance independently. This will mean
            % incorporating the wt fxns seperately for the mean and variance.
            
            % (1) make a spatiotemporal wt fxn. The template will peak at
            % +/- 1 instead of the peak Rstar.
            switch params.enableScones
                case true
                    wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
                case false
                    wt_spaceTime = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, 0);
            end
            
            % filter the wt fxn through the cone IRF
            [wt_spaceTime, ~] = coneVolution_FFT(wt_spaceTime, gab, cones, [0 0 0]);
            
            % the function that convolves the wtFxn also adjusted the output
            % according to the cone gain. This is a hack to undo this
            % normalization so that the wt fxn peaks/troughs at +/-1.for
            % each cone type.
            scaleFact = max(max(max(abs(wt_spaceTime),[],1),[],2),[],3);
            wt_spaceTime = bsxfun(@rdivide, wt_spaceTime, scaleFact);         
            
            % New code added 8/9/19 GDLH to change temporal integration
            t = [0:1./cones.samplingRate:(gab.framesInNewGabor-1)/cones.samplingRate];
            zeroedTimeSteps = t>max(params.temporalIntegrationLims) | t<min(params.temporalIntegrationLims);
            wt_spaceTime(:,:,zeroedTimeSteps,:) = zeros(size(wt_spaceTime(:,:,zeroedTimeSteps,:)));
            idlob.wt_spaceTime = wt_spaceTime;
            % End of new code added 8/9/19 GDLH
            
            % (2) calculate the mean response of the ideal observer for a
            % prototypical cone in each pixel location. Since the model
            % averages signals across cones with in the same pixel, the
            % 'average' response is identical to the 1 cone prototypical
            % response.
            idlOb_mean = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            idlOb_mean = idlOb_mean .* wt_spaceTime;
            idlOb_mean = sum(idlOb_mean,3); % summing across time is part of the 'dot product' with the weight function
            
            % (2.5) Sum signals across retinal locations (i.e., pixels).
            % Exclude pixels that are aperatured. Normally, none of the
            % pixels are aperatured, but for some comparisons with V1,
            % there are a few pixels missing.
            cone_mask{1} = cones.num_L > 0;
            cone_mask{2} = cones.num_M > 0;
            cone_mask{3} = cones.num_S > 0; 
            for a = 1:3
                tmp = squeeze(idlOb_mean(:,:,:,a)); % X by Y for a single cone type (the 'time' dim was summed across above).;
                idlob.analyticMean{clr, cnt}(1,a) = reshape(tmp, 1, []) * cone_mask{a}(:); % a row vector times a column vector.
            end

            
            % (3) calculate the variance of the ideal observer.
            wt_spaceTime_FFT = fft(wt_spaceTime,[],3);
            noise_fft = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_fft, [1,3,2])); % weighting by the space/time/clr fxn in the fft domain
            
            % normalizing by N due to discrete fourier transform (see
            % Plancherel's theorem)
            N = size(wt_spaceTime, 3);
            noise_fft = noise_fft./N;
            
            % finish the 'dot product' in the fourier domain by summing
            % across frequencies
            noise_fft = sum(abs(noise_fft).^2, 3); % variance of a single pixel across freqs
            noise_fft = squeeze(noise_fft); % remove the (singelton) freq dimension
            
            % pool across cones (within a pixel) and then across
            % pixels. Exclude pixels due to the aperature (when
            % present). The ideal observer sums signals from adjecent
            % pix, so variances add. This is accomplished in the code
            % by summing across pixel locations. Lastly, the noise is
            % independent from cone to cone (within pixels) so the
            % average var = var(x)/numCones
            numCones{1} = cones.num_L;
            numCones{2} = cones.num_M;
            numCones{3} = cones.num_S;
            for a = 1:3
                tmp_noise = reshape(noise_fft(:,:,a), [], 1);
                tmp_num = numCones{a}(:);
                idx = cone_mask{a}(:);
                idlob.analyticVar{clr, cnt}(1,a) = sum(tmp_noise(idx) ./ tmp_num(idx)); % noise is independent, so just divide by the number of cones.
            end
            
        case 'equate_mosaic'
            
            % For cone isolating stimuli, the ideal observer should produce
            % the identical responses if the cone mosaics are ths same. Use
            % this functionality to run some debugging code.
            cones.num_M = cones.num_L;
            cones.num_S = cones.num_L;
    end
end

function [idlob, cones] = obsMethod_phaseInvariant(command, idlob, gab, cones, mon, clr, cnt, trl, params)
    
    switch lower(command)
        case 'initialize'
            nColors = size(gab.colorDirs,1);
            nContrasts = size(gab.contrasts{1},2);
            nTrials = gab.nTrials;

            % dimensionality will be (nColors x nContrasts x nTrials) for
            % the monte carlo method, and (nColors x nContrasts) for the
            % analytic solution. Each entry will be a CELL ARRAY of six
            % numbers, one for each cone type and each wt fxn.
            [idlob.resp(1:nColors, 1:nContrasts, 1:nTrials)] = deal({nan(1,6)});
            [idlob.analyticMean(1:nColors, 1:nContrasts), idlob.analyticVar(1:nColors, 1:nContrasts)] = deal({nan(1,6)});


        case 'compute'
            
            % subtract off the mean from the linear response
            zeroMeanLinResp = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            
            % generate the pooled response and add the (pooled) noise.
            sigPlusNoise = zeroMeanLinResp + pooledNoise_average(cones.noise, cones.num_L, cones.num_M, cones.num_S);
            clear zeroMeanLinResp %kill to save ram
            
            phases = [0, pi/2]; % sin and cos
            for a = 1:2 %loop over the two wt fxns.
                
                % weight the sigPlusNoise by a template based on the gabor
                % stimulus. The template will peak at +/- 1 instead
                % of the peak Rstar. The template cone mosaics are all IN PHASE
                % despite the fact that the stimulus may not be.
                switch params.enableScones
                    case true
                        wt = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, phases(a));
                    case false
                        wt = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, phases(a));
                end
                [wt, ~] = upsampleAndBookend(wt, cones, [0 0 0]);
                tRange = (cones.nFrontPad+1) : (size(wt,3) - cones.nBackPad);
                wt = wt(:,:, tRange, :); %hack off the bookends.
                
                dotWithLinResp = sum(sigPlusNoise .* wt, 3); %sum across the 'time' dimension
                idealObsResp = sum(sum(dotWithLinResp,1),2); %sum across the spatial dimensions.
                idx = [((a-1).*3)+1 : a.*3];
                idlob.resp{clr, cnt, trl}(idx) = squeeze(idealObsResp);
            end
            
        case 'analytic_solution'
            % compute the mean and variance independently. This will mean
            % incorporating the wt fxns seperately for the mean and variance.
            phases = [0 pi/2];
            for a = 1:2;
                
                % index for storing things in the idlob array
                idx = [((a-1).*3)+1 : a.*3];
                
                % (1) make a spatiotemporal wt fxn. The template will peak at
                % +/- 1 instead of the peak Rstar.
                switch params.enableScones
                    case true
                        wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, phases(a));
                    case false
                        wt_spaceTime = makeGaborMovie([1,1,0], gab, mon, 'idlob template', cones, params, phases(a));
                end
                [wt_spaceTime, ~] = upsampleAndBookend(wt_spaceTime, cones, [0 0 0]);
                tRange = (cones.nFrontPad+1) : (size(wt_spaceTime,3) - cones.nBackPad);
                wt_spaceTime = wt_spaceTime(:,:, tRange, :); %hack off the bookends.
                
                
                % (2) calculate the mean response of the ideal observer
                idlOb_mean = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
                idlOb_mean = idlOb_mean .* wt_spaceTime;
                idlOb_mean = sum(sum(sum(idlOb_mean,1),2),3); %sum across space. The sum over time is actually part of the dot product with the wt. fxn...
                idlob.analyticMean{clr, cnt}(idx) = squeeze(idlOb_mean);
                
                % (3) calculate the variance of the ideal observer. This
                % should be the same for each contrast (within color
                % directions) so just calculate it on the first trial and then
                % carry it forward for the other contrasts
                if cnt == 1;
                    
                    % only do the fft once.
                    if a == 1;
                        wt_spaceTime_FFT = fft(wt_spaceTime,[],3);
                    end
                    
                    % apply the wt fxn to the noise PS.
                    noise_fft = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_fft, [1,3,2])); % weighting by the space/time/clr fxn in the fft domain
                    
                    % normalizing by N due to discrete fourier transform (see
                    % Plancherel's theorem)
                    N = size(wt_spaceTime, 3);
                    noise_fft = noise_fft./N;
                    
                    % finish the 'dot product' in the fourier domain by summing
                    % across frequencies
                    noise_fft = sum(abs(noise_fft).^2, 3); % variance of a single pixel across time
                    noise_fft = squeeze(noise_fft); % remove the (singelton) time dimension
                    
                    % pool across cones (within a pixel) and then across pixels
                    numCones = cat(3, cones.num_L, cones.num_M, cones.num_S);
                    noise_fft = noise_fft ./ numCones; % the noise is independent from cone to cone (within pixels) so the average var = var(x)/numCones
                    noise_fft = sum(sum(noise_fft,1),2); % ideal observer sums signals from adjecent pix, so variances add.
                    
                    idlob.analyticVar{clr, cnt}(idx) = squeeze(noise_fft);
                    
                else
                    idlob.analyticVar{clr, cnt}(idx) = idlob.analyticVar{clr, 1}(idx);
                end
            end
            
        case 'equate_mosaic'

            % For cone isolating stimuli, the ideal observer should produce
            % the identical responses if the cone mosaics are ths same. Use
            % this functionality to run some debugging code.
            cones.num_M = cones.num_L;
            cones.num_S = cones.num_L;
    end
    
end

function [idlob, cones] = obsMethod_absThresh(command, idlob, gab, cones, mon, clr, cnt, trl, params)
    %
    % This ideal observer only has access to the spatiotemporal weighting
    % function. Also, the space/time wt function is a gabor stimulus that
    % peaks/troughs at �1. All the cones in the wt fxn are in phase with
    % eachother. The result is a 3 vector for each trial. I'll have to do
    % linear discriminant inorder to determine the neurometric threshold.
    
    switch lower(command)    
        case 'initialize'
            nColors = size(gab.colorDirs,1);
            nContrasts = size(gab.contrasts{1},2);
            nTrials = gab.nTrials;
            
            % dimensionality will be (nColors x nContrasts x nTrials) for
            % the monte carlo method, and (nColors x nContrasts) for the
            % analytic solution. 
            [idlob.resp] = deal(nan(nColors, nContrasts, nTrials));
            [idlob.analyticMean, idlob.analyticVar] = deal(nan(nColors, nContrasts));
            
        case 'compute'
            % subtract off the mean from the linear response
            zeroMeanLinResp = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            
            % generate the pooled response and add the (pooled) noise.
            sigPlusNoise = zeroMeanLinResp + pooledNoise_average(cones.noise, cones.num_L, cones.num_M, cones.num_S);
            clear zeroMeanLinResp %kill to save ram
            
            % weight the sigPlusNoise by a template based on the cone's
            % temporal impulse response. Make sure that the wt fxn peaks at
            % one.
            IRF = ifft(cones.filter_fft);
            IRF = IRF ./ max(IRF);
            IRF = IRF(1:cones.filterSamples);
            wt = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0); % use [1,1,1] to make the SAME template for all color directions in the simulation
            [wt, ~] = upsampleAndBookend(wt, cones, [0 0 0]);
            tRange = (cones.nFrontPad+1) : (size(wt,3) - cones.nBackPad);
            wt = wt(:,:, tRange, :); %hack off the bookends.
            wt(:) = 0; % for absThresh, turn all the values to zero.
            wt(1,1,:,1) = permute(IRF, [2, 4, 1, 3]);  
            
            dotWithLinResp = sum(sigPlusNoise .* wt, 3); %sum across the 'time' dimension
            idealObsResp = sum(sum(dotWithLinResp,1),2); %sum across the spatial dimensions.
            idlob.resp(clr, cnt, trl) = idealObsResp(1,1,1,1);
            
        case 'analytic_solution'
            % compute the mean and variance independently. This will mean
            % incorporating the wt fxns seperately for the mean and variance.
            
            % (1) make a spatiotemporal wt fxn. The template will peak at
            % +/- 1 instead of the peak Rstar.
            IRF = ifft(cones.filter_fft);
            IRF = IRF ./ max(IRF);
            IRF = IRF(1:cones.filterSamples);
            wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
            [wt_spaceTime, ~] = upsampleAndBookend(wt_spaceTime, cones, [0 0 0]);
            tRange = (cones.nFrontPad+1) : (size(wt_spaceTime,3) - cones.nBackPad);
            wt_spaceTime = wt_spaceTime(:,:, tRange, :); %hack off the bookends.
            wt_spaceTime(:) = 0;
            wt_spaceTime(1,1,:,1) = permute(IRF, [2, 4, 1, 3]);
            
            % (2) calculate the mean response of the ideal observer
            idlOb_mean = bsxfun(@minus, cones.linresp, permute(mon.bkgndlms_pA, [3, 2, 4, 1])); %make the linear response zero mean.
            idlOb_mean = idlOb_mean .* wt_spaceTime;
            idlOb_mean = sum(sum(sum(idlOb_mean,1),2),3); %sum across space. The sum over time is actually part of the dot product with the wt. fxn...
            idlob.analyticMean(clr, cnt) = idlOb_mean(1,1,1,1); % since there's only cones in the Lmosaic, just take the L cone data
            
            % (3) calculate the variance of the ideal observer.
            wt_spaceTime_FFT = fft(wt_spaceTime,[],3);
            noise_fft = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_fft, [1,3,2])); % weighting by the space/time/clr fxn in the fft domain
            
            % normalizing by N due to discrete fourier transform (see
            % Plancherel's theorem)
            N = size(wt_spaceTime, 3);
            noise_fft = noise_fft./N;
            
            % finish the 'dot product' in the fourier domain by summing
            % across frequencies
            noise_fft = sum(abs(noise_fft).^2, 3); % variance of a single pixel across time
            noise_fft = squeeze(noise_fft); % remove the (singelton) time dimension
            
            % pool across cones (within a pixel) and then across pixels
            numCones = cat(3, cones.num_L, cones.num_M, cones.num_S);
            noise_fft = noise_fft ./ numCones; % the noise is independent from cone to cone (within pixels) so the average var = var(x)/numCones
            idlob.analyticVar(clr, cnt) = noise_fft(1,1,1,1); % since there's only cones (in one pix) in the Lmosaic, only take data from this pix and do not sum signals across the other pixels
    end
end

function pooledNoise = pooledNoise_average(noise, num_L, num_M, num_S, DEBUGMODE)
    
    %determine the number of each cone type per pixel    
    numCones = nan([size(num_L), 1, 3]);
    numCones(:,:,1,1) = num_L;
    numCones(:,:,1,2) = num_M;
    numCones(:,:,1,3) = num_S;    

    % In debug mode, store the original variances
    if (nargin==5) && strcmpi(DEBUGMODE, 'unit_test')
        old_var = var(noise, [], 3);
        pool_prediction = bsxfun(@rdivide, old_var, numCones);
    end


    % Calculate the pooled signal and noise. For i.i.d draws:
    % Var(mean(X,Y,...Z)) = N*Var ./ N^2
    % Var(a .* coneNoise) = a^2 .* Var(coneNoise)
    % If a^2 * Var(x) = N*Var ./ N^2.
    % Than  a = sqrt(1/N).
    % For i.i.d draws, the mean of means can be approximated by the sample mean 
    pooledNoise = bsxfun(@times, sqrt(1./numCones), noise); % calculate the variance after pooling
    
    % In debug mode, show that the noise after pooling is what you'd expect
    % based on the pre-pooling variance
    if exist('DEBUGMODE', 'var') && strcmpi(DEBUGMODE, 'unit_test')
        new_var = var(pooledNoise, [], 3);
        figure, hold on,
        pltClr = {'r', 'g', 'b'};
        sym = {'o', '.', 's'};
        for a = 1:3
            plot(pool_prediction(:,:,a), new_var(:,:,a), [sym{a}, pltClr{a}], 'markersize', 10)
        end
        XX = get(gca, 'xlim');
        plot([0, XX(2)], [0, XX(2)], 'k')
        title('Pooled Var: predited vs. computed')
        xlabel('Predicted Variance')
        ylabel('Actual Variance')
        drawnow
    end
end

function [cones, mon] = setupAbsThreshParams(cones, mon, params)
    % This function sets things up for the absolute threshold experiments.
    
    % We want to restrict our attention to only a few cones, so take the
    % original cone mosaic and modify it accordingly.
    Ncones = params.Ncones; % the TOTAL number of cones to consider.
    [cones.num_L, cones.num_M, cones.num_S] = deal(zeros(size(cones.num_L))); % start by zeroing everything out.
    cones.num_L(1,1) = Ncones; % put all the cones in the same pix. To make life easy, just put it in 1st pix.
    
    % Change the bkgnd adaptation state of the cones. Re-calculate the cone
    % gain scale factor and the bkgnd photo current
    mon.bkgndlms_Rstar = [0;0;0];
    cones = defineConeGainScaleFactor(cones, mon);
    mon = getBkgndLinearResp(cones, mon);
end



function prefDiam = getAperatureTuning(stro, spikenum)
    
    % CAH: I stole this code from getGratingsTuning.m
    
    diams = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'diam'));
    protocols = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'protocol'));
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_on'));
    stimoff_t= stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_off'));


    spikerates = [];
    baseline_t = 0.25;
    for i = 1:size(stro.trial,1)
        spiketimes = stro.ras{i,spikenum};
        nspikes = sum(spiketimes > stimon_t(i) & spiketimes < stimoff_t(i));
        spikerates = [spikerates; nspikes./(stimoff_t(i)-stimon_t(i))];
        nspikes = sum(spiketimes > stimon_t(i)-baseline_t & spiketimes < stimon_t(i));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Looking at area summation curve
    if (any(protocols == 5))
        Lprotocol = protocols == 5;
        x = diams(Lprotocol);
        y = spikerates(Lprotocol);
        pp = csape(x,y,'not-a-knot');
        xx = linspace(min(diams),max(diams),100);
        fit = ppval(pp,xx);
        prefDiam = xx(find(fit == max(fit),1));

    else
        prefDiam = nan;
    end



end



function saveDataToDisk(gab, cones, mon, idlob, params) %#ok<INUSL>
    fprintf('  Saving the data...\n')
    
    % grab the contents of the .m file used to generate the data (i.e.,
    % this file you're reading now)
    pathToDTcones = which('DTcones');
    fid = fopen(pathToDTcones, 'r');
    mfiletext = fscanf(fid, '%c');
    fclose(fid);
    
    %cd to where the data should be stored
    originalDir = pwd;
    newDir = params.saveDir;
    if isempty(newDir)
        newDir = fullfile(strtok(userpath, ':'), 'DTcones');
    end
    cd(newDir)
    
    % make a directory. Make a new one if this is a simple call to DTcones,
    % or make a specific directory (which will be deleted later) if this is
    % called during a parfor loop with a dtnt run.
    if params.parallelOperations
        switch lower(params.runType)
            case 'dtnt'
                sepIdx = regexp(params.DTNT_fname, filesep);
                timestamp = params.DTNT_fname(sepIdx(end)+1 : end);
            case 'dtv1'
                [~, timestamp, ~] = fileparts(params.DTV1_fname);
                
        end
    else
        d = dir;
        suffix = 1;
        timestamp = [date, sprintf('_%d', suffix)];
        alreadyExists = any(strcmpi(timestamp, {d.name}));
        while alreadyExists
            suffix = suffix+1;
            timestamp = [date, sprintf('_%d', suffix)];
            alreadyExists = any(strcmpi(timestamp, {d.name}));
        end
    end
    
    %save the data
    cones.sigPlusNoise = []; % delete the arrays that are unnecessary (and large)
    cones.linresp = [];
    gab.movie = [];
    mkdir(pwd, timestamp)
    cd(timestamp)
    fname_dat = ['out_', timestamp, '.mat'];
    eval(['save ' fname_dat, ' gab cones mon idlob params'])
    
    %save the code of the calling function (i.e., DTcones.m)
    fname_mfile = ['DTcones.m_', timestamp];
    fid = fopen(fname_mfile, 'w');
    fwrite(fid, mfiletext, 'char');
    fclose(fid);
    
    %be nice and cd back to where you came from
    cd(originalDir)
end

function out = softEq(A, B, precision, dim)

%
%   SOFT EQUALS
%
%  EXAMPLE: softEquals = softEq(A, B, precision, dim);
%
% Determines the element-wise equality between matricies A and B but
% neglects any round off error. The amount error to be tolerated is
% specified by the optional 'precision' argument. The default tolerance is
% ten decimal places. The optional argument `dim` specifies whether the
% equality should be computed along the 'rows' or 'cols'.
%
% CAH 03/09

if nargin < 2 || isempty(A) || isempty(B)
    out = [];
    return
end

if nargin < 3 || isempty(precision);
    precision = 10;
end

%if the optional argument dim is specified then treat this as a call to
%ismember
if nargin > 3
    if strcmpi(dim, 'cols')
        A = A';
        B = B';
    end
    
    %make sure A is the smaller of the two matricies. I'll use A as the
    %template and look into B for instances of A.
    if size(A,1) > 1
        error(['the template (first argument) must be a row or column vector' ...
            ' (depending on what you passed in as `dim`)'])
    end
    matches = abs(B - repmat(A, size(B, 1), 1)) < 10^(-precision);
    out = sum(matches, 2) == size(B, 2);
    
    %transform back if the caller specified columns
    if strcmpi(dim, 'cols')
        out = out';
    end
else
    %just the element wise soft equals:
    out = abs(A-B) < 10^(-precision);
end
end

%
%   FUNCTIONS USED FOR TESTING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function callAllTestingSubfunctions(gab, mon, cones, params, idlob)

    % redefine the default colors of the gabor (for testing purposes)
    gab.colorDirs = [1 1 1; 1 -1 0];
    gab.contrasts = [{.17}, {.10}];
    
    %
    % run the following testing routines.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    test_compareRstarAndRGBgabors(gab, mon, params, cones);
    test_GaborMovie(gab, mon, params, cones)
    test_coneMosaic(cones);
    test_coneVolution_FFTvsConv(gab, mon, cones, params);
    test_shapedConeNoise(cones, gab, mon, params);
    test_impulseResponse(cones, gab, mon, params);
    test_fftWeightFunction(cones, gab, mon, params)
    test_phaseInvariance(gab, cones, mon, params)
    
    fprintf('\n\n\n********** TESTING COMPLETE ************\n')
end

function test_GaborMovie(gab, mon, params, cones)

    % This testing unit just plots the Gabor movie as a sanity check for
    % the spatial frequency, orientation, and overall R*/sec rates for
    % stimulus and bkgnd.

    fprintf('\n\n\n *********** Gabor Movie Test *************\n')
    fprintf('  Theta: %.3f radians\n', gab.theta)
    fprintf('  Spt Freq: %.3f cpd\n', gab.sf)
    fprintf('  Sigma: %.3f tenths DVA\n', gab.sd)
    
    clr = 1;
    cnt = 1;
    LMS_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    fprintf('  Color dir = [%s]\n', num2str(gab.colorDirs(clr,:)));
    
    figure
    set(gcf, 'position', [21 430 1411 392]);
    fontsize = 16;
    fontname = 'helvetica';
    
    
    for b = 1:2:size(gab.movie,3)
        subplot(1,3,1)
        tmp = permute(gab.movie(:,:,b,:), [1,2,4,3]);
        L = tmp(:,:,1);
        imagesc(L)
        axis square
        colorbar
        set(gca, 'fontsize', fontsize, 'fontname', fontname);
        title('R* L-Cone Image')
        
        subplot(1,3,2)
        M = tmp(:,:,2);
        imagesc(M)
        axis square
        colorbar
        set(gca, 'fontsize', fontsize, 'fontname', fontname);
        title('R* M-Cone Image')
        xlabel(sprintf('FRAME NUMBER: %d', b))
        
        subplot(1,3,3)
        S = tmp(:,:,3);
        imagesc(S)
        axis square
        colorbar
        set(gca, 'fontsize', fontsize, 'fontname', fontname);
        title('R* S-Cone Image')
        
        %update as fast as possible
        drawnow
    end
    
    
    %show a single frame of Lcc - Mcc image to see how they add. For L-M
    %stimuli, the non-zero component of the sum should be very small
    clr = 2;
    cnt = 1;
    LMS_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    
    figure
    frame = round(size(gab.movie, 3)./2);
    tmp = permute(gab.movie(:,:,frame,:), [1,2,4,3]);
    Lcc = (tmp(:,:,1) - mon.bkgndlms_Rstar(1))./mon.bkgndlms_Rstar(1);
    Mcc = (tmp(:,:,2) - mon.bkgndlms_Rstar(2))./mon.bkgndlms_Rstar(2);
    imagesc(Lcc+Mcc); %in CC space
    axis square
    colorbar
    set(gca, 'fontsize', fontsize, 'fontname', fontname);
    title('Lcc map + Mcc map (should be low amp for L-M stim)')
    xlabel('Pix X')
    ylabel('Pix Y')
    drawnow
end

function test_compareRstarAndRGBgabors(gab, mon, params, cones)
    %I got nervous that computing the R* for a trial b/4 computing the
    %gabor movie would yield different results than computing the movie in
    %rgb space and then tranforming to R* space. This function compares the
    %two techniques
    
    %determin the rgb/R* to test
    rgb_example = mon.bkgndrgb + [0.4, -0.4, 0];
    LMS_Rstar = mon.rgb2Rstar * rgb_example(:);
    
    %generate the gabors
    gab_Rstar = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    gab_rgb = makeGaborMovie_rgb(rgb_example, gab, mon);
    
    figure
    hold on,
    maxVal = max([gab_rgb.movie(:); gab_Rstar(:)]);
    minVal = min([gab_rgb.movie(:); gab_Rstar(:)]);
    plot(gab_rgb.movie(:), gab_Rstar(:), '.')
    plot([minVal-100, maxVal+100], [minVal-100, maxVal+100], 'k-')
    title('Gabor_rgb vs. Gabor_Rstar')
    xlabel('Gabor prepared with rgb then transformed to R*')
    ylabel('Gabor prepared directly from R*')
    drawnow
    
    fprintf('\n\n\n********** comparing rgb and R* gabors ************\n')
    fprintf('  All rgb gabor equals R* gabor: %d\n', all(abs(gab_rgb.movie(:) - gab_Rstar(:))<1e-12))
    dif = LMS_Rstar - mon.bkgndlms_Rstar;
    fprintf('  LMS R* = [%.3f, %.3f, %.3f] (diff from bkgnd)\n', dif(1), dif(2), dif(3))
    
    %create a local version of the gabor movie function, but make it operate on
    %rgbs, and then convert to R* (which is opposite the way the main code
    %works)
    function gab = makeGaborMovie_rgb(rgb, gab, mon)

        % Define the constants
        nPix = round((gab.nSd*2) .* (gab.sd./10) .* mon.pixperdeg);
        halfSize = round(nPix./2);
        gaborPhase = 0;
        flashDeltaPhase = gab.driftRate * 2 * pi * (1./mon.frameRate); %the amount to advance each frame
        sigmaInPix = (gab.sd./10) .* mon.pixperdeg;
        lambdaInPix = 1./(gab.sf./mon.pixperdeg);
        rgb_increment = rgb - mon.bkgndrgb;

        % Parameters for ramping on/off stimulus

%         flashNumFrames = ceil(mon.frameRate * gab.length);
%         rampLength = ceil(flashNumFrames / 4); %in frames
%         ramp = linspace(0, 1, rampLength);  %ramp is 1/4th of the total duration on either side
%         plateau = ones(1,flashNumFrames - (rampLength .* 2));
%         flashTimeProf = [ramp, plateau, fliplr(ramp)];

        % meshgrid for the gabor
        row = -halfSize:halfSize-1; %subtract one so that you don't overlap the texture window;
        col = -halfSize:halfSize-1;
        [X, Y] = meshgrid(row, col);
        xprime = X .* cos(-gab.theta) + Y .* sin(-gab.theta);
        yprime = -X .* sin(-gab.theta) + Y .* cos(-gab.theta);


        % Preallocate space for the movie (4D: [X, Y, Time, Guns])
        gab.movie = nan(size(X,1), size(X,1), flashNumFrames, 3);

        for t = 1:flashNumFrames;

            % Multiply the gabor by the increment from background and by the temporal weighting fxn.
            gabor = exp(-(xprime.^2 + gab.gamma.^2 .* yprime.^2) ./ (2.* sigmaInPix.^2)) .* cos(2 .* pi .* yprime ./ lambdaInPix + gaborPhase);

            %gun intensities b/w 0 and 1
            gab.movie(:,:,t,1) = (gabor .* gab.flashTimeProf(t) .* rgb_increment(1)) + mon.bkgndrgb(1);
            gab.movie(:,:,t,2) = (gabor .* gab.flashTimeProf(t) .* rgb_increment(2)) + mon.bkgndrgb(2);
            gab.movie(:,:,t,3) = (gabor .* gab.flashTimeProf(t) .* rgb_increment(3)) + mon.bkgndrgb(3);

            % Update the phase
            gaborPhase = gaborPhase + flashDeltaPhase;
        end

        %now convert to R* using some complicated reshaping of the gabor matrix
        tmp = reshape(gab.movie, [], 3)';
        tmp = mon.rgb2Rstar * tmp; %converts from rgb to R*
        gab.movie = reshape(tmp', size(gab.movie,1), size(gab.movie,2), flashNumFrames, 3);
    end
end

function test_coneMosaic(cones)

    % Test to check the distribution of cones across the mosaic

    fprintf('\n\n\n********** plotting the cone mosaic ************\n')
    fprintf('  L/M cones should be equally likely\n')
    fprintf('  S-cone prevalance should be ~10 percent \n')

    N = cones.num_L + cones.num_M + cones.num_S;
    
    figure
    set(gcf, 'position', [8 305 1424 524])
    colormap(pmkmp(250, 'CubicL'))
    fontsize = 14;
    
    subplot(2,3,1) %l cone map
    imagesc(cones.num_L)
    colorbar
    set(gca, 'fontsize', fontsize)
    title('Number of L cones by pixel')
    xlabel('pix X')
    ylabel('pix Y')
    
    subplot(2,3,4) % histo of num_L cones per pix
    hist(cones.num_L(:)./N(:))
    set(gca, 'fontsize', fontsize)
    xlabel('Percent L cones per pixel')
    ylabel('Counts')
    
    subplot(2,3,2) %m cone map
    imagesc(cones.num_M)
    colorbar
    set(gca, 'fontsize', fontsize)
    title('Number of M cones by pixel')
    xlabel('pix X')
    ylabel('pix Y')
    
    subplot(2,3,5) % histo of num_M cones per pix
    hist(cones.num_M(:)./N(:))
    set(gca, 'fontsize', fontsize)
    xlabel('Percent M cones per pixel')
    ylabel('Counts')
    
    subplot(2,3,3) %s cone map
    imagesc(cones.num_S)
    colorbar
    set(gca, 'fontsize', fontsize)
    title('Number of S cones by pixel')
    xlabel('pix X')
    ylabel('pix Y')
    
    subplot(2,3,6) % histo of num_S cones per pix
    hist(cones.num_S(:)./N(:))
    set(gca, 'fontsize', fontsize)
    xlabel('Percent S cones per pixel')
    ylabel('Counts')
    
    drawnow
end

function test_coneVolution_FFTvsConv(gab, mon, cones, params)

    % A test that the FFT-based convolution is the same as the time-based
    % convolution
    
    fprintf('\n\n\n********** comparing FFT-convolution to the real convolution ************\n')
    fprintf('  This should take a moment... \n')
    
    % should the stimulus be an impulse or a regular gabor?
    STIMULUSTYPE = 'gabor'; %for impulse response, set to 'impulse'
    fprintf('  The stimulus type is: %s\n', STIMULUSTYPE);
    
    % make a gabor movie from the default colors at high contrast:
    clr = 1;
    cnt = 1;
    LMS_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
    fprintf('  Color dir = [%s]\n', num2str(gab.colorDirs(clr,:)));
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    [cones.linresp, gab.movie] = upsampleAndBookend(gab.movie, cones, mon.bkgndlms_Rstar);
    tmpMovie = cones.linresp; %put the gabor movie into a temporary array
    
    if strcmpi(STIMULUSTYPE, 'impulse')
        tmpMovie(:,:,:,:) = 0;
        idx = ((cones.convLength-cones.nFrontPad-cones.nBackPad)/2) + cones.nFrontPad;
        tmpMovie(:,:,idx,:) = cones.samplingRate;
    end
    
    % resize the movie and filter to the dimesnions appropriate for
    % convolution
    filter = ifft(cones.filter_fft);
    filter = filter(1:(cones.nBackPad+1)); %should be the appropriate length
    filter = filter(:);
    inds = (cones.convLength-cones.nBackPad+1): cones.convLength;
    tmpMovie(:,:,inds,:) = [];
    
    % do the convolution for a random subset of the pixels/cones
    pixIdx = unique(unidrnd(size(tmpMovie,1), [16,2]), 'rows');
    coneIdx = unidrnd(3, size(pixIdx,1), 1);
    resp_conv = {};
    for a = 1:size(pixIdx,1) %do the comparison for some randomly selected pixels and cones
        timevec = tmpMovie(pixIdx(a,1), pixIdx(a,2), :, coneIdx(a));
        timevec = timevec(:);
        tmp = conv(filter, timevec)./cones.samplingRate;
        tmp = tmp .* cones.gainScaleFactor(coneIdx(a));
        resp_conv{a} = tmp(cones.nFrontPad+1 : (cones.convLength-cones.nBackPad));
    end
    
    % now do the FFT convolution and compare the resuts to the regular
    % convolution
    cones.linresp = [];
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    
    if strcmpi(STIMULUSTYPE, 'impulse')
        [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar, 'impulse_response');
    else
        [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar);
    end
    
    resp_fft = {};
    for a = 1:size(pixIdx,1)
        tmp = cones.linresp(pixIdx(a,1), pixIdx(a,2), :, coneIdx(a));
        resp_fft{a} = tmp(:);
    end
    
    %plot the results
    nPlts = ceil(sqrt(size(pixIdx,1)));
    figure
    set(gcf, 'name', 'CONV (black) and FFT (red)')
    for a = 1:size(pixIdx,1)
        subplot(nPlts, nPlts, a)
        hold on,
        plot(resp_fft{a}, 'ro')
        plot(resp_conv{a}, 'k', 'linewidth', 2)
        axis tight
        title(sprintf('Pix: [%d, %d], Cone: %d', pixIdx(a,1), pixIdx(a,2), coneIdx(a)))
    end
    drawnow
    
    %plot the results (diff b/w time and fft conv)
    nPlts = ceil(sqrt(size(pixIdx,1)));
    figure
    set(gcf, 'name', 'CONV - FFT')
    for a = 1:size(pixIdx,1)
        subplot(nPlts, nPlts, a)
        hold on,
        plot(resp_conv{a}-resp_fft{a}, 'b')
        title(sprintf('Pix: [%d, %d], Cone: %d', pixIdx(a,1), pixIdx(a,2), coneIdx(a)))
        axis tight
        hold off
    end
    drawnow
end

function test_shapedConeNoise(cones, gab, mon, params)
    
    % Testing various aspects of the cone-noise generation.
    
    fprintf('\n\n\n********** checking computations of cone noise ************\n')
    fprintf(' This should take a moment... \n')
    
    % make a gabor movie from the default colors at high contrast:
    fprintf('  Setting up a test Gabor stimulus... \n')
    clr = 1;
    cnt = 1;
    LMS_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
    fprintf('  Color dir = [%s]\n', num2str(gab.colorDirs(clr,:)));
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);
    
    % calculate the linear response of the cones
    fprintf('  Calculating the linear response... \n')
    [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar);
    
    % now call shapedConeNoise with the optional "testing" argument
    fprintf('  Now adding cone noise... \n')
    cones = shapedConeNoise(cones, gab, 'unit_test');
    
    
    % see if the average variance (in time) is the same as the integral of
    % the power spectrum
    N = size(cones.noise,3);
    var_tt = var(cones.noise,[],3);
    var_ff = sum(cones.modelNoise_ps)./(N.*(N-1));
    [counts, bins] = hist(var_tt(:), 20);
    figure, hold on,
    bar(bins, counts)
    plot([var_ff, var_ff], [0, max(counts)], 'r', 'linewidth', 2)
    drawnow
    
    %demonstrate that caculating the pooled variance gives what you would
    %expect baeed on the variance of the time series before pooling.
    fprintf('  Caculating the pooled variance... \n')
    sigPlusNoise = cones.linresp + pooledNoise_average(cones.noise, cones.num_L, cones.num_M, cones.num_S, 'unit_test');
    
end

function test_impulseResponse(cones, gab, mon, params)

    fprintf('\n\n\n*********** Impulse Response Function Test *************\n')
    fprintf('   compare the model''s impulse response to the simulated IRF\n')
    
    % Simulating the impulse response
    %
    % calculate the linear response to a unit step in Rstar. First,
    % generate a random movie, and then pass it to coneVolution with the
    % optional 'unit_test' argument. This will make a discrete delta
    % function in Rstar from a zero bkgnd. The linear response to the delta
    % function should be identical to Juan's empirically derived impulse
    % response.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gab.movie = makeGaborMovie(mon.bkgndlms_Rstar, gab, mon, 'movie', cones, params, 0);
    [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar, 'impulse_response');
    
    %plot the impulse response    
    figure, hold on,
    stimOnIdx = round(gab.framesInNewGabor./2);
    coneExamp = permute(cones.linresp(1,1, :, 1), [3,1,2,4]);
    coneExamp = coneExamp(stimOnIdx:end);
    tt = ((1:numel(coneExamp))-1)./cones.samplingRate;
    irf = ifft(cones.filter_fft); % I don't think I need to normalize this b/c the filter_fft wasn't normalized...
    irf = irf(1:numel(coneExamp));
    irf = irf .* cones.gainScaleFactor(1); % assuming that I'm plotting the L-cone response
    plot(tt, irf, 'ro')
    plot(tt, coneExamp, 'b', 'linewidth', 2)
    set(gca, 'fontsize', 14)
    xlabel('time (sec)')
    ylabel('pico Amps per R star')
    title('Temporal impulse response')
    legend('simulated', 'model')
    drawnow
    
    % TWO: 40 R* in 10ms on 4000 R*/sec Background
    %
    % This is an experiment that juan performed on real cones, and provides
    % a point of comparisons between the actual cones and the model cones.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('   next, test the response to a 4000 R*/sec pulse on a 4000 R*/sec bkgnd\n')
    
    PULSEWIDTH = 0.010; % seconds
    PULSEDELAY = 0.200; % seconds after stim on to position the pulse
    PULSEBKGND = 4000;  % in R*/sec
    PULSEHEIGHT = 40;   % in R*... NOT R*/sec!!!!
    
    
    % convert the pulseheight to R*/sec during the pulse... and to
    % accomodate for the new sampling rate.
    effectiveRStarPerSec = ((PULSEBKGND .* PULSEWIDTH) +  PULSEHEIGHT) ./ PULSEWIDTH
    effectivePulseWidth = ceil(PULSEWIDTH.*mon.frameRate) .* (cones.samplingRate./mon.frameRate) ./ cones.samplingRate
    
    % make a gabor movie and then change it to make the pulse step.
    gab.movie = makeGaborMovie(mon.bkgndlms_Rstar, gab, mon, 'movie', cones, params, 0);
    gab.movie(:) = PULSEBKGND;
    tvec = [0:size(gab.movie,3)-1] ./ mon.frameRate;
    ton_idx = find(tvec > PULSEDELAY,1, 'first');
    nFramesOn = ceil(PULSEWIDTH .* mon.frameRate);
    gab.movie(:,:,[ton_idx : (ton_idx+nFramesOn-1)],:) = effectiveRStarPerSec;
    
    % re-estimate the cone linear response to the new bkgnd.
    mon.bkgndlms_Rstar(:) = PULSEBKGND;
    cones = defineConeGainScaleFactor(cones, mon);
    cones.gainScaleFactor
    mon = getBkgndLinearResp(cones, mon);
    
    % run the conevolutions and estimate the cone noise
    [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar);
    cones = shapedConeNoise(cones, gab);
    resp = cones.linresp + cones.noise;   
    
    
    % plot some of the results
    nExamps = 10; % the number of example traces
    pixels = [unidrnd(size(cones.noise,1), nExamps,1), unidrnd(size(cones.noise,2), nExamps,1), unidrnd(size(cones.noise,4), nExamps,1)];
    space = 20;
    cone_tvec = [0:(size(cones.linresp,3)-1)] ./ cones.samplingRate;
    
    % single trial examples
    figure, hold on
    for a = 1:nExamps
        tmp = permute(resp(pixels(a,1), pixels(a,2), :, pixels(a,3)), [3,1,2,4]);
        yy = ((a-1).*space) + tmp;
        plot(cone_tvec, yy, 'k')
    end
    % add the example stimulus
    t_pulse = [ton_idx : (ton_idx+nFramesOn)]./ mon.frameRate;
    yy = (mean(tmp)+(a*space)) .* ones(1,numel(t_pulse));
    plot(t_pulse, yy, 'b', 'linewidth', 3)
    xlabel('Time (sec)')
    ylim([mean(tmp)-space, yy(1)+10])
    xlim([0, (size(resp,3)-1) ./ cones.samplingRate])
    set(gcf, 'position', [397 5 447 824])
    %set(gca, 'yticklabel', [])
    title(sprintf('Response to %d R* pulse', PULSEHEIGHT))
    drawnow

    % plot the mean of the cone signal+noise
    figure, hold on,
    set(gca, 'fontsize', 14)
    r = [];
    for a = 1:3
        tmp = resp(:,:,:,a);
        tmp = permute(tmp, [3, 1, 2]);
        tmp = reshape(tmp, size(resp,3), []);
        r = [r,tmp];
    end
    r = mean(r,2);
    plot(cone_tvec, r, 'k');
    yy = (max(r).*1.01) .* ones(1,numel(t_pulse));
    plot(t_pulse, yy, 'b', 'linewidth', 3)
    xlim([0, cone_tvec(end)])
    xlabel('Time (sec)')
    ylabel('Pico Amps')
    drawnow
    
end

function test_fftWeightFunction(cones, gab, mon, params)

    % This testing module concerns the application of the spaceTime weight
    % function to the noise array during the analytical solution of the
    % pooled response. This script tries to trouble shoot
    % applying the wt fxn in the Fourier domain.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %
    % First a toy problem: The PS of a noise process should be equal to the
    % diagonal elements of the COV(fft coeffs).
    %
    %%%%%%%%
    sigma = 4.5;
    N = 1000;
    noise = normrnd(0, sigma, 2000, N);
    
    % cov in freq domain
    noise_fft = fft(noise, [], 2);
    cov_ff = cov(noise_fft);
    ps_from_cov = fftshift(diag(cov_ff));
    
    % pow spect (Avg)
    noise_ps = abs(noise_fft.^2);
    ps_from_fft = fftshift(mean(noise_ps, 1));
    
    % variance of fft coeff
    var_from_fft = fftshift(var(noise_fft, [], 1));
    
    figure, hold on,
    plot(ps_from_fft, 'b')
    plot(ps_from_cov, 'bo')
    plot(var_from_fft, 'r+')
    
    
    %
    % Set things up to test the actual portions of the model. make a linear
    % response and noise array.
    %
    %%%%%%%
    
    fprintf('\n\n\n*********** Shaped Noise Test *************\n')
    fprintf('   compare different implementations of shaping the cone noise spectra and deriving the variance \n')
    
    
    clr = 1; cnt = 1;
    fprintf('  Setting up a test Gabor stimulus... \n')
    fprintf('  Color dir = [%s]\n', num2str(gab.colorDirs(clr,:)));
    LMS_Rstar = getGaborRstar(clr, cnt, gab, mon, params, cones);
    gab.movie = makeGaborMovie(LMS_Rstar, gab, mon, 'movie', cones, params, 0);

    fprintf('  Calculating the linear response... \n')
    [cones.linresp, gab.movie] = coneVolution_FFT(gab.movie, gab, cones, mon.bkgndlms_Rstar);

    fprintf('  Now adding cone noise... \n')
    cones = shapedConeNoise(cones, gab);
    
    %
    % make the weight function that is the same for every pixel. 
    %
    %%%%%%%%
    wt_spaceTime = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
    [wt_spaceTime, ~] = upsampleAndBookend(wt_spaceTime, cones, mon.bkgndlms_Rstar);
    tRange = (cones.nFrontPad+1) : (size(wt_spaceTime,3) - cones.nBackPad);
    wt_spaceTime = wt_spaceTime(:,:, tRange, :); %hack off the bookends.
    
    dummy = ones(size(wt_spaceTime));
    nX = round(size(dummy,1)/2);
    wt_spaceTime = bsxfun(@times, dummy, wt_spaceTime(nX, nX, :, 1));
    Nsamps = size(wt_spaceTime, 3);    
    
    %
    % plot a histogram of variances for each time series, and that
    % predicted by the model's power spectrum
    %
    %%%%%%%%%%%
    N = size(cones.linresp,3);
    normfact_ff = 1./(N.*(N-1));
    var_preWeight_tt = var(cones.noise, 0, 3);
    var_preWeight_ff = sum(cones.modelNoise_ps).*normfact_ff;
    figure, hold on,
    hist(var_preWeight_tt(:), 20)
    plot([var_preWeight_ff, var_preWeight_ff], [0, 100], 'r', 'linewidth', 3)
    
    
    %
    % "weight" the noise array in the time domain.
    %
    %%%%%%%%
    noise_timewt_tt = cones.noise .* wt_spaceTime;  
    

    
    %
    % "weight" the noise in the Fourier domain using a simple piece-wise
    % mulitiplication
    %
    %%%%%%%%
    wt_spaceTime_FFT = fft(wt_spaceTime,[],3);
    noise_fftwt_ff = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_fft, [1,3,2])); % weighting by the space/time/clr fxn in the fft domain
    
    noise_pswt_ff = bsxfun(@times, wt_spaceTime_FFT, permute(cones.modelNoise_ps, [1,3,2]));
    noise_pswt_ff = noise_pswt_ff .* wt_spaceTime_FFT;
    var_from_ps = unique(sum(abs(noise_pswt_ff),3));
    var_from_ps = var_from_ps ./ N^2 %Plancherel's theorm says that [x_t' * y_t] => [1/N (X_f' * Y_f)] but here i'm using a PS, so normalize by N^2?
    
    %
    % Estimate the covariance matrix for the time and frequency
    % representations of the shaped noise
    %
    %
    %%%%%%%%%%%%%%%%
    shapedNoise = permute(cones.noise, [3,1,2,4]); % time goes down columns
    shapedNoise = reshape(shapedNoise, Nsamps, [])'; % time goes across rows
    cov_noise_tt = cov(shapedNoise);
    
    shapedNoise_ff = fft(shapedNoise,[],2);
    cov_noise_ff = cov(shapedNoise_ff);
    
    figure
    subplot(1,2,1)
    imagesc(cov_noise_tt)
    colorbar
    
    subplot(1,2,2)
    imagesc(abs(cov_noise_ff))
    colorbar
    
    figure, hold on,
    plot(cones.freqAxis_fft, diag(cov_noise_ff), 'k.')
    plot(cones.freqAxis_fft, cones.modelNoise_ps, 'b')
    
    
    %
    % Now the part that actually matters: The analytic solution needs to
    % solve for the variance of the dot-products in the time domain.
    %
    %%%%%%%%%
    dotproducts = sum(noise_timewt_tt,3);
    dotproducts = dotproducts(:);
    var_dotprod = var(dotproducts)
    
    % prediction from the fft method. Divide by N b/c Plancherel's theorem
    % say's you have to...
    var_fftwt_ff = unique(sum(abs(noise_fftwt_ff./N).^2, 3))
   
    
    
    figure, hold on,
    bins = linspace(min(dotproducts), max(dotproducts), 50);
    counts = histc(dotproducts, bins);
    counts = counts./sum(counts);
    bar(bins, counts, 'type', 'histc')
    pred_old = normpdf(bins, mean(dotproducts), sqrt(var_fftwt_ff));
    plot(bins, pred_old./sum(pred_old), 'r', 'linewidth', 2)
    drawnow
    
    
    % finish with another toy problem: just demonstrate that the method
    % that we use to generate the variance of the dot product is legit.N = 100;
    
    % two random arrays (a wt fxn and a "noise ps")
    iters = 1000;
    N = 100;
    wt = normrnd(0,1,iters,N);
    noise = normrnd(0,1,iters,N);
    
    dot_in_time = diag(wt*noise');
    
    dot_in_fft = diag((1/N) .* fft(wt,[],2)*ctranspose(fft(noise,[],2)));
    
    figure %demonstrate that the time version and the fft version are the same
    plot(dot_in_time, dot_in_fft, '.')
    
    
    % now try from the PS of x, which is the way I do it in the model
    noise_ps = abs(fft(noise,[],2)).^2;
    
    tmp = fft(wt,[],2) .* sqrt(noise_ps);
    tmp = abs(tmp).^2;
    var_of_dot = (1/N^2) .* sum(tmp, 2);
    
    mean(var_of_dot)
    var(dot_in_time)
   
end

function test_phaseInvariance(gab, cones, mon, params)


    fprintf('\n\n\n *********** Test Phase Invariant (6D) idlob method *************\n')

    
    
    oldTF = gab.driftRate;
    plotTF = 3; % more plots for this specific TF.
    %tf = [0, plotTF, logspace(log10(0.00001), log10(40), 40)];
    tf = [0, 0.01, 0.5, plotTF, 20, 40];
    tf = sort(tf);

    [xbar_raw, xbar_proj1, xbar_proj2] = deal(nan(size(tf)));
    for i = 1:numel(tf);
        
        fprintf('TF %d of %d \n', i, numel(tf));
        
        % reassign the TF
        gab.driftRate = tf(i);
        
        % make the sin and cos wt fxns.
        wt_sin = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, 0);
        [wt_sin, ~] = upsampleAndBookend(wt_sin, cones, [0 0 0]);
        tRange = (cones.nFrontPad+1) : (size(wt_sin,3) - cones.nBackPad);
        wt_sin = wt_sin(:,:, tRange, :); %hack off the bookends.
        wt_sin = wt_sin(:); % turn it into a column vector.
        wt_sin = wt_sin ./ norm(wt_sin); % make it a unit vector.
        
        wt_cos = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, pi/2);
        [wt_cos, ~] = upsampleAndBookend(wt_cos, cones, [0 0 0]);
        wt_cos = wt_cos(:,:, tRange, :); %hack off the bookends.
        wt_cos = wt_cos(:);
        wt_cos = wt_cos ./ norm(wt_cos);

        % make a bunch of phase shifted versions of the wt fxn, for each
        % version, project out the sin and cos wt fxn, and keep track of the
        % result after projecting it out...
        niters = 25;
        phases = linspace(0, 2*pi, niters);

        [beforeProj, projOut_cosAndSin, projOut_cos] = deal(nan(niters, 1));
        for a = 1:niters
            
            % make the phase shifted version
            wt_raw = makeGaborMovie([1,1,1], gab, mon, 'idlob template', cones, params, phases(a));
            [wt_raw, ~] = upsampleAndBookend(wt_raw, cones, [0 0 0]);
            wt_raw = wt_raw(:,:, tRange, :); %hack off the bookends.
            wt_raw = wt_raw(:);
            beforeProj(a) = norm(wt_raw);

            % project out the wt_cos
            projVal = wt_raw' * wt_cos;
            wt_raw = wt_raw - (projVal .* wt_cos);
            projOut_cos(a) = norm(wt_raw);

            % project out the wt_sin from the result
            projVal = wt_raw' * wt_sin;
            wt_raw = wt_raw - (projVal .* wt_sin);
            projOut_cosAndSin(a) = norm(wt_raw);

        end
        
        xbar_raw(i) = median(beforeProj);
        xbar_proj1(i) = median(projOut_cos);
        xbar_proj2(i) = median(projOut_cosAndSin);
        
        
        if tf(i) == plotTF;
            figure
            subplot(1,3,1)
            hist(beforeProj);
            title('raw wt fxn (all phases)')
            ylabel('Counts')
            
            subplot(1,3,2)
            hist(projOut_cos);
            title('proj out wt\_cos')
            xlabel('Norm of result')
            
            subplot(1,3,3)
            hist(projOut_cosAndSin);
            title('proj out wt\_sin + wt\_cos')
        end

    end
    
    figure, hold on,
    idx = tf == 0;
    pltTF = tf;
    pltTF(idx) = min(pltTF(~idx)) .* .2;
    
    
    plot(pltTF, xbar_raw, 'k.-')
    plot(pltTF, xbar_proj1, 'b.-')
    plot(pltTF, xbar_proj2, 'r.-')
    set(gca, 'xscale', 'log', 'yscale', 'log')
    legend('raw', 'proj cos', 'proj both')
    
    
    % reassign the TF
    gab.driftRate = oldTF;
    
end






