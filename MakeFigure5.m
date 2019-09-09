% Script for generating Figure 5

load ('LGN_data_stro');

MMPERDEG = 0.223; % mm/deg (Perry and Cowey 1985)
DEGPERMM = 1/MMPERDEG; % deg/mm
MONKEYS={'Monkey 1','Monkey 2'};
CELLTYPES = {'M','P'};
TEMPORONASALSCALEFACTOR = .8;
ONOFFCORRELATION = 0.05;
RFTRUNCATIONINSD = 2;
HUMAN2MONKPSCALEFACTOR = .80; % From Dacey and Petersen. reasonable range: [.77 .81];
POPULATIONSCALING = 1;

ecc_to_diam_deg_M = @(rf_r_deg) 10.^(-1.2459+0.0345*rf_r_deg); % temporal retina equivalent
a = 0.9729; % Table 1
r2 = 1.084; % Table 1
re = 7.633; % Table 1
dc_0 = 14804.6; % Cone density of fovea
rm = 41.03; % See Equation 7
ecc_to_diam_deg_P = @(x)(sqrt(2./(sqrt(3).*... % Equation 9. Distance between adjacent midget RF centers.
    (2*dc_0.*(1+x./rm).^-1.*(a*(1+(x./r2)).^-2+(1-a)*exp(-x./re)))...
    ./2))... % Dividing y by 2 to estimate RF size from only ON or OFF mosaics (halving the density).
    *HUMAN2MONKPSCALEFACTOR); % Monkey midget RFs are slightly smaller than human midget RFs

% Stuff for DTcones_gh
fundamentals = load('T_cones_smj10');
params = [];
params.runType = 'isosamp';
params.obsMethod = 'obsMethod_filteredWtFxn';
params.impulseResponse = 'rieke';
params.DTV1_fname = [];
params.DTNT_fname = [];
params.unitTest = false;
params.eqMosaic = false;
%params.saveDir = '/Users/greghorwitz/Documents/MATLAB';
%params.notes = 'IsoSamp test';
params.parallelOperations = false;
params.eyeType = 'monkey';
params.eyeNumber = 1; % LGN neurons are monocular
params.coneSampRate = 2400;
params.flatPowerSpect = false;
params.enableScones = false;

out = cell(2,2);
for monkey_idx = 1:size(data,1)
    for celltype_idx = 1:size(data,2)
        for i = 1:length(data{monkey_idx,celltype_idx})
            stro = data{monkey_idx, celltype_idx}{i};
            sigmas_n = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigmas_n')));
            gabor_sigma = unique(stro.trial(:,strcmp(stro.sum.trialFields(1,:),'sigma')));

            [tmp.uniquestim, tmp.dprime] = IsoSampGetDPrime(stro);
            
            Lblank = tmp.uniquestim(:,1) == 0 & tmp.uniquestim(:,2) == 0 & tmp.uniquestim(:,3) == 0;
            Llum = sign(tmp.uniquestim(:,1)) == sign(tmp.uniquestim(:,2)) & ~Lblank;
            
            % Now getting population scale factor
            % Getting RF size at rfx, rfy
            if strcmp(CELLTYPES{celltype_idx},'M')
                ecc_to_diam_deg = ecc_to_diam_deg_M;
            else
                ecc_to_diam_deg = ecc_to_diam_deg_P;
            end
            if POPULATIONSCALING
                [tmp.population_scalefactor, tmp.ncells] = IsoSampGetPopulationScaleFactor(stro, ecc_to_diam_deg, TEMPORONASALSCALEFACTOR,RFTRUNCATIONINSD, ONOFFCORRELATION, 2)
            else
                tmp.population_scalefactor = 1; tmp.ncells = 1;
            end
            rfx = stro.sum.exptParams.rf_x/10;
            rfy = stro.sum.exptParams.rf_y/10;
            tmp.ecc = sqrt((rfx./TEMPORONASALSCALEFACTOR)^2+rfy^2);

            % Now getting cone ideal observer sensitivity
            % Cut and paste from IsoSampPop section 4.1
            params.stro = stro;
            spds = params.stro.sum.exptParams.mon_spd;
            spds = SplineSpd([380:4:780]',reshape(spds,length(spds)/3,3),[380:5:780]');
            M = fundamentals.T_cones_smj10*spds;
            cal.monSpect = spds(:);
            cal.Mmtx = M(:);
            cal.frameRate = params.stro.sum.exptParams.framerate;
            cal.bkgndrgb = params.stro.sum.exptParams.bkgndrgb';
            cal.fname = 'test';
            cal.monSpectWavelengths = linspace(380,780,101);
            cal.pixperdeg = params.stro.sum.exptParams.pixperdeg;
            params.monCalFile = cal;
            if POPULATIONSCALING
                params.gab.sd = gabor_sigma;
                params.eyeNumber = 2;
            else
                params.gab.sd = ecc_to_diam_deg(tmp.ecc)/2; % overwriting gabor SD with RF SD
                params.gab.sd = params.gab.sd; % DEBUGGING: artificially inflating RF SD by 5x
                params.eyeNumber = 1;
            end
            
            [gab, cones, mon, idlob, params] = IsoSampConeCurrents(params,0);
            
            conedata = [];
            for k = 1:size(idlob.analyticMean,1) % looping over color direction
                for j = 1:size(idlob.analyticMean(k,:),2) % looping over contrast/TF
                    if ~isempty(idlob.analyticMean{k,j})
                        tmp_lm_mu = idlob.analyticMean{k,j}([1 2]);
                        tmp_lm_var = idlob.analyticVar{k,j}([1 2]);
                        tf = gab.driftRates{k}(j);
                        conedata = [conedata; gab.colorDirs(k,[1 2]).*gab.contrasts{k}(j) tf tmp_lm_mu tmp_lm_var];
                    end
                end
            end
            
            % Using uniquestim to order the rows of tmpdata and calculating cone
            % dprimes
            tmp.conedprime = [];
            for j = 1:size(tmp.uniquestim,1)
                L = all(abs(tmp.uniquestim(j,:)-conedata(:,[1 2 3]))<1e-10,2);
                if sum(L) ~= 1
                    if all(tmp.uniquestim(j,:) == 0)
                        tmp.conedprime(j) = nan;
                    else
                        error('sorting problem');
                    end
                else
                    v = conedata(L,[6 7]); % variance
                    m = conedata(L,[4 5]); % mean
                    tmp.conedprime(j) = sqrt(m.^2*(1./v'));
                end
            end
            
            % Need to calculate flashTimeProfile on the basis of timewindoffset
            t = 0:1/mon.frameRate:gab.length;
            nframes = length(t);
            flashTimeProfile = ones(1,nframes);
            ramp = linspace(0,1,nframes/4);
            flashTimeProfile(1:length(ramp)) = ramp;
            flashTimeProfile(end:-1:end-length(ramp)+1) = ramp;
            sigmaInPix = params.gab.sd*cal.pixperdeg;
            
            tmp.photondprime = IsoSampGetPhotonDPrime(flashTimeProfile, mon.frameRate, mon.bkgndlms_Rstar, sigmaInPix, cat(3,cones.num_L,cones.num_M), tmp.uniquestim);
            % Adding to the "data" cell array of cell arrays
            listsofar = out{monkey_idx, celltype_idx};
            listsofar{length(listsofar)+1} = tmp;
            out{monkey_idx, celltype_idx}=listsofar;
        end
    end
end

% Defining a few constants for needed for the TF binning
TFBINCENTERS = logspace(log10(1), log10(60), 15);
binwidth = log10(TFBINCENTERS(2))-log10(TFBINCENTERS(1));
TFBINEDGES = 10.^(linspace(log10(TFBINCENTERS(1))-binwidth/2,log10(TFBINCENTERS(end))+binwidth/2,length(TFBINCENTERS)+1));

% Doing the plotting
figure('Position',[440 100 750 700],'DefaultAxesTickDirMode','manual','DefaultAxesTickdir','out','DefaultAxesYcolor','black','DefaultAxesXcolor','black')
set(gcf,'DefaultAxesFontSize',15,'DefaultAxesFontAngle','italic','DefaultAxesUnits','centimeters');
hax = [];
hax(1,1) = axes('position',[2 2 9 9]); hold on;
hax(1,2) = axes('position',[2 14 9 9]); hold on;
hax(2,1) = axes('position',[14 2 9 9]); hold on;
hax(2,2) = axes('position',[14 14 9 9]); hold on;

for monkey_idx = 1:size(data,1)
    for celltype_idx = 1:size(data,2)
        monkey_celltype_data = out{monkey_idx,celltype_idx};
        tmp = nan*ones(4,length(TFBINCENTERS),length(monkey_celltype_data)); % rows: LGN dprime, cone dprime, photon dprime, ns
        for i = 1:length(monkey_celltype_data)
            cell_data = monkey_celltype_data{i};
            Lblank = cell_data.uniquestim(:,1) == 0 & cell_data.uniquestim(:,2) == 0 & cell_data.uniquestim(:,3) == 0;
            Llum = sign(cell_data.uniquestim(:,1)) == sign(cell_data.uniquestim(:,2)) & ~Lblank;
            for j = 1:length(TFBINCENTERS)
                Ltf = cell_data.uniquestim(:,3) > TFBINEDGES(j) & cell_data.uniquestim(:,3) <= TFBINEDGES(j+1);
                if sum(Ltf&Llum) > 0
                    tmp(1,j,i) = mean(cell_data.dprime(Llum&Ltf))*cell_data.population_scalefactor;
                    tmp(2,j,i) = mean(cell_data.conedprime(Llum&Ltf));
                    tmp(3,j,i) = mean(cell_data.photondprime(Llum&Ltf));
                    tmp(4,j,i) = 1; % each cell counts as an independent entity
                end
            end
        end
        allnans = all(isnan(squeeze(tmp(4,:,:))),2);
        tmp = tmp(:,~allnans,:);
        n = nansum(tmp(4,:,:),3);
        mn_lgn = nanmean(tmp(1,:,:),3);
        sd_lgn = nanstd(tmp(1,:,:),0,3);
        sem_lgn = sd_lgn./sqrt(n);
        mn_cone = nanmean(tmp(2,:,:),3);
        sd_cone = nanstd(tmp(2,:,:),0,3);
        sem_cone = sd_cone./sqrt(n);
        mn_photon = nanmean(tmp(3,:,:),3);
        sd_photon = nanstd(tmp(3,:,:),0,3);
        sem_photon = sd_photon./sqrt(n);
        
        axes(hax(monkey_idx,celltype_idx));
        TBC = TFBINCENTERS(1:length(n));
        patch([TBC, fliplr(TBC)],[mn_photon+sem_photon, fliplr(mn_photon-sem_photon)],[0 .7 .7],'Facealpha',.5);
        plot(TBC, mn_photon,'c-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','cyan');
        patch([TBC, fliplr(TBC)],[mn_cone+sem_cone, fliplr(mn_cone-sem_cone)],[1 .5 .5],'Facealpha',.5);
        plot(TBC, mn_cone,'r-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','red');
        patch([TBC, fliplr(TBC)],[mn_lgn+sem_lgn, fliplr(mn_lgn-sem_lgn)],[.5 .5 .5],'Facealpha',.5);
        plot(TBC, mn_lgn,'k-o','LineWidth',1,'MarkerEdgeColor','black','MarkerFaceColor','black');
        plot([TFBINEDGES(1) TFBINEDGES(end)],[1.27 1.27],'k-');
        if POPULATIONSCALING
            set(gca,'Xscale','log','Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'Ylim',[-2 40]);
        else
            set(gca,'Xscale','log','Xlim',[TFBINEDGES(1) TFBINEDGES(end)],'Ylim',[-2 10]);
        end
        xlabel('Frequency (Hz)');
        ylabel('SNR');
        title([MONKEYS{monkey_idx},' ',CELLTYPES{celltype_idx}])
    end
end
set(gcf,'Renderer','painters');