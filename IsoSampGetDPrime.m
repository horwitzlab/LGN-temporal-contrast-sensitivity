function [uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro, method, whichspike, signaloffsets, TESTING)
    % [uniquestim, dprime, signal, noise] = IsoSampGetDPrime(stro, method, [whichspike], [signaloffsets], [TESTING]);
    %
    % INPUTS
    %      stro: an IsoSamp stro structure.
    %      method: how to calculate d'
    %           1) Pooling variance of noise and signal+noise within TF condition (default)
    %           2) Pooling noise variance across TFs, not using
    %           signal+noise variance at all (not a great choice)
    %           3) Pooling noise variance and signal+noise variance across TF conditions.
    %
    % Takes an IsoSamp stro structure as an input and returns d-prime values
    % for each stimulus condition.
    %
    % Optional arguments: whichspike, signaloffsets, TESTING
    %     whichspike = 1 or 2 (for sig001a and sig001b, respectively)
    %     signaloffsets: (1x2) matrix. Offset in seconds for spike counting
    %     relative to stimon and stimoff respectively. [0 0] is the
    %     default.
    %     TESTING just makes an additional plot
    %
    % GDLH 6/18/17
        
    if nargin < 2
        method = 1;
    end
    if ~ismember(method,[1 2 3])
        error('Unknown noise distribution estimate method (1 2 or 3)');
    end
    if nargin < 3
       spikestring = 'sig001a';
    elseif whichspike == 2
        spikestring = 'sig001b';
    else
        spikestring = 'sig001a';
    end
    
    stimon_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimon_t'));
    stimoff_t = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stimoff_t'));

    if nargin < 4 || isempty(signaloffsets)
        signaloffsets = [0 mean(stimoff_t-stimon_t)];
    end
    dur = signaloffsets(2)-signaloffsets(1);
    if nargin < 5
       TESTING = 0;
    end
    if stro.sum.paradigmID ~= 107
        error('Input argument not an IsoSamp stro structure');
    end
        
    Lcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_l'));
    Mcc = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'stim_m'));
    TF = stro.trial(:,strcmp(stro.sum.trialFields(1,:),'tf'));
    spikeidx = strcmp(stro.sum.rasterCells(1,:),spikestring);
    spikes = stro.ras(:,spikeidx);
    uniquestim = sortrows(unique([Lcc Mcc TF],'rows'),3); % sorting by TF
    if nargout == 1
        return
    end
    
    totaltrials = size(stro.trial,1); 
    spiketimes = {}; % First, getting all the spike times
    nspikes = [];
    Lblank = Lcc==0 & Mcc==0 & TF==0;
    if sum(Lblank) == 0
        disp(['No blank trials in file ',stro.sum.fileName]);
    end
    signal_data = nan(totaltrials,2); % dot products onto basis vectors at stimulus fundamental frequency
    
    % Getting dot products
    for i = 1:totaltrials
        current_TF = TF(i);
        tmpspikes = spikes{i}-stimon_t(i);
        spiketimes{i} = tmpspikes(tmpspikes > signaloffsets(1) & tmpspikes < signaloffsets(end));
        spiketimes{i} = spiketimes{i} - signaloffsets(1);
        nspikes(i) = length(spiketimes{i});
        signal_data(i,:) = [sum(cos(2*pi*current_TF.*spiketimes{i})) sum(sin(2*pi*current_TF.*spiketimes{i}))];
    end
    
    % Calculating Mahalanobis distances, first for stimulus present trials.
    % For details of the calculation, see section 7 of StatsStuff.m
    T2s = [];
    for i = 1:totaltrials
        T2s(i) = getStats(signal_data(i,:),nspikes(i),TF(i),dur);
    end
    % Getting an empirical distribution of Mahalanobis distances (per TF)
    % for the zero contrast trials. This is in a separate loop because 
    % we're using the
    % same noise trials for each TF.

    empirical_noise_T2 = [];
    for i = 1:size(uniquestim,1)
        current_TF = uniquestim(i,3);
        tmp = [];
        for j = find(Lblank)'
            noiseprojections = [sum(cos(2*pi*current_TF.*spiketimes{j})) sum(sin(2*pi*current_TF.*spiketimes{j}))];
            tmp = [tmp; getStats(noiseprojections,nspikes(j),current_TF,dur)];
        end
        empirical_noise_T2{i} = tmp;
    end
    
    % When there's only one spike, it's always the same distance away from
    % the origin = no variance in T2.
    
    % "Signal" and "noise" distributions as a function of condition
    signal = {};
    noise = {};
    for i = 1:size(uniquestim,1)
        L = Lcc == uniquestim(i,1) & Mcc == uniquestim(i,2) & TF == uniquestim(i,3);
        signal{i} = norminv(chi2cdf(T2s(L),2),0,1);
        signal{i}(isinf(signal{i})) = nan; % T2s = 0 converted to snr = nan. Ignored by nanmean, etc.
        noise{i} = norminv(chi2cdf(empirical_noise_T2{i},2),0,1);
        noise{i}(isinf(noise{i})) = nan; % T2s = 0 converted to snr = nan. Ignored by nanmean, etc.
    end
    
    % Debugging
    if TESTING
        tmp = [];
        figure; subplot(2,1,1); hold on;
        for i = 2:length(signal) % Skipping the zero contrast condition
            errorbar(uniquestim(i,3),nanmean(signal{i}),nanstd(signal{i}),'k-');
            plot(uniquestim(i,3),nanmean(signal{i}),'ko','markerfacecolor','black');
            errorbar(uniquestim(i,3),nanmean(noise{i}),nanstd(noise{i}),'b-');
            plot(uniquestim(i,3),nanmean(noise{i}),'bo','markerfacecolor','blue');
            tmp = [tmp; noise{i}];
            %            plot(repmat(uniquestim(i,3),size(signal{i})),signal{i},'ko');
        end
        set(gca,'Xscale','log','xlim',[.5 100]);
        title(stro.sum.fileName(find(stro.sum.fileName == '/',1,'last')+1:end));
        subplot(2,1,2);
        hist(tmp,[-4:.1:4]);
    end
    dprime = zeros(length(signal),1);
    if method == 1 % pooling noise variance and signal+noise variance within TF condition
        for i = 1:length(signal)
            n_signal = sum(~isnan(signal{i}));
            n_noise = sum(~isnan(noise{i}));
            pooled_var = ((n_signal-1)*nanvar(signal{i})+(n_noise-1)*nanvar(noise{i}))/(n_signal+n_noise-2);
            dprime(i) = (nanmean(signal{i})-nanmean(noise{i}))/sqrt(pooled_var);
        end
    elseif method == 2 % Pooling noise variance across TFs, not using signal+noise variance
        tmp = [];
        for i = 1:length(noise)
            tmp = [tmp; noise{i}];
        end
        noisemean = nanmean(tmp);
        noisevar = nanvar(tmp);
        for i = 1:length(signal)
            dprime(i) = (nanmean(signal{i})-noisemean)/sqrt(noisevar);
        end
    else % pooling noise variance and signal+noise variance across TF conditions. Might want to use pooled variance estimate.
        tmp = [];
        for i = 1:length(noise)
            tmp = [tmp; nanvar(noise{i}); nanvar(signal{i})];
        end
        for i = 1:length(signal)
            dprime(i) = (nanmean(signal{i})-nanmean(noise{i}))/sqrt(nanmean(tmp));
        end
    end
    
    % Nested function that returns the mean and variance of the projections
    % of Poisson spike trains onto truncated sine and cosine functions and
    % the Mahalanobis distance to the empirical spike train.
    function [dM2, mu, S2] = getStats(projs,nspikes,TF,dur)
        mu = nspikes*[sin(2*pi*TF*dur)+0 -cos(2*pi*TF*dur)+1]/(2*pi*TF*dur);
        
        term1 = (nspikes/8)*(sin(4*pi*TF*dur)/(pi*TF*dur)+4);
        term2 = (nspikes*(nspikes-1)/dur^2)*sin(2*pi*TF*dur)^2/(2*pi*TF)^2;
        crossprods(1,1) = term1+term2-mu(1)^2;
        
        term1 = (nspikes*sin(2*pi*TF*dur)^2/(4*pi*TF*dur));
        term2 = (nspikes*(nspikes-1)/dur^2)*sin(pi*TF*dur)^3*cos(pi*TF*dur)/(pi*TF)^2;
        crossprods(2,1) = term1+term2-prod(mu);
        crossprods(1,2) = crossprods(2,1);
        
        term1 = (nspikes/dur)*((dur/2)-(sin(4*pi*TF*dur)/(8*pi*TF)));
        term2 = (nspikes*(nspikes-1)/dur^2)*sin(pi*TF*dur)^4/(pi*TF)^2;
        crossprods(2,2) = term1+term2-mu(2)^2;
        S2 = crossprods;
        if isnan(rcond(S2)) || rcond(S2) == 0
            dM2 = 0;
        else
            dM2 = (projs-mu)*inv(S2)*(projs-mu)';
        end
    end
end