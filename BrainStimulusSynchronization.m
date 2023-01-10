clear
close all

fs = 512;
sub = 63:81;
sub([5,15]) = [];

load('Chan128.mat'); % channel locations
chanlocs = cell(1,128);
for i=1:length(chanlocs)
    chanlocs{i} = chan(i).labels;
end

%% 
spec = cell(length(sub),2);
SigS001 = cell(length(sub),1);
SigS002 = cell(length(sub),1);
for s=1:length(sub)   
    address = ['\N',num2str(sub(s)),'\PreFFT.mat']; % Load preprocessed data of each subject
    load(address);
    
% There are two types of events (Duple/Triple Rhythm (S001) and Quadruple
% Rhythm (S002)
    LS001 = 0;
    LS002 = 0;
    for i=1:length(EEG.event)
        if strcmp(EEG.event(i).type,'S001')
            LS001 = LS001 + 1;
        elseif strcmp(EEG.event(i).type,'S002')
            LS002 = LS002 + 1;
        end
    end
    
    sigS001 = zeros(EEG.nbchan,floor(38*EEG.srate),LS001);
    sigS002 = zeros(EEG.nbchan,floor(36*EEG.srate),LS002);
    k1 = 0;
    k2 = 0;
    for i=1:length(EEG.event)
        if strcmp(EEG.event(i).type,'S001')
            k1 = k1 + 1;
            sigS001(:,:,k1) = EEG.data(:,EEG.event(i).latency:...
                EEG.event(i).latency + floor(38*EEG.srate)-1);
        elseif strcmp(EEG.event(i).type,'S002')
            k2 = k2 + 1;
            sigS002(:,:,k2) = EEG.data(:,EEG.event(i).latency:...
                EEG.event(i).latency + floor(36*EEG.srate)-1);
        end
    end
    
    tmp = cell(length(EEG.chanlocs),1);
    for i=1:length(EEG.chanlocs)
        tmp{i} = EEG.chanlocs(i).labels;
    end
    temp = ismember(chanlocs,tmp);
    tmp1 = find(temp==0);
    tmp2 = find(temp==1);
    NewSig(tmp2,:,:) = sigS001;
    NewSig(tmp1,:,:) = NaN;
    sigS001 = NewSig;    
    clear NewSig
    NewSig(tmp2,:,:) = sigS002;
    NewSig(tmp1,:,:) = NaN;
    sigS002 = NewSig;    
    clear NewSig
    
    SigS001{s} = sigS001;
    SigS002{s} = sigS002;                
end

%% Wavelet transform on EEG data
cfg = [];
cfg.channel    = 'all';
cfg.method     = 'wavelet';
cfg.pad        = 'nextpow2';
cfg.width      = 7;
cfg.output     = 'pow';
cfg.foi        = .5:.25:10;
cfg.gwidth     = 3;
cfg.polyremoval = 0;
cfg.padtype = 'zero';

time = linspace(0,38,38*EEG.srate);
for s=1:length(sub)
    data = SigS001{s};
    data(:,:,sum(isnan(squeeze(sum(data,2))))==128) = [];
    S001 = zeros(length(chan),length(cfg.foi),38*EEG.srate);
    for i=1:size(data,3)
        cfg.toi = linspace(10,48,38*EEG.srate);
        temp = data(:,1:length(time),i);      
        
        dat = [zeros(size(temp,1),10*fs),temp,zeros(size(temp,1),10*fs)]; 
        padding = 2^nextpow2(length(dat));
        cfg.pad = padding/fs;
        options = {'pad', cfg.pad, 'padtype', cfg.padtype, 'freqoi', cfg.foi,...
            'polyorder', cfg.polyremoval};
        [spectrum,foi,toi] = ft_specest_wavelet(dat, linspace(0,58,length(dat)),...
            'timeoi', cfg.toi, 'width', cfg.width, 'gwidth', cfg.gwidth,options{:});
        spectrum(:,:,1:2*fs) = [];
        spectrum(:,:,end-2*fs+1:end) = [];

        S001 = S001 + spectrum./abs(spectrum);
    end  
    spec{s,1} = S001./size(data,3);
end

time = linspace(0,36,36*EEG.srate);
for s=1:length(sub)
    data = SigS002{s};
    data(:,:,sum(isnan(squeeze(sum(data,2))))==128) = [];
    S002 = zeros(length(chan),length(cfg.foi),36*EEG.srate);
    for i=1:size(data,3)
        cfg.toi = linspace(10,46,36*EEG.srate);
        temp = data(:,1:length(time),i);      
        
        dat = [zeros(size(temp,1),10*fs),temp,zeros(size(temp,1),10*fs)]; 
        padding = 2^nextpow2(length(dat));
        cfg.pad = padding/fs;
        options = {'pad', cfg.pad, 'padtype', cfg.padtype, 'freqoi', cfg.foi, 'polyorder', cfg.polyremoval};
        [spectrum,foi,toi] = ft_specest_wavelet(dat, linspace(0,56,length(dat)), 'timeoi', cfg.toi, 'width', cfg.width, 'gwidth', cfg.gwidth,options{:});
        spectrum(:,:,1:4*fs) = [];
        spectrum(:,:,end-4*fs+1:end) = [];
       
        S002 = S002 + spectrum./abs(spectrum);
    end  
    spec{s,2} = S002./size(data,3);
end

%% Wavelet transform on Stimuli
fs = 512;
Wn = 100;

[RhytI,fT]  = audioread('\Duple-Triple Rhythm.wav');
RhytI       = RhytI(:,1);
RhytI       = envelope(RhytI,Wn,'analytic');
r1          = 10;
RhytI       = decimate(RhytI,r1,'fir'); % downsample from  44100 Hz to 4410 Hz
r2          = round(4410/512);
RhytI       = decimate(RhytI,r2,'fir'); % downsample from 4410Hz to 490Hz
p           = 256;
q           = 245;
RhytI       = resample(RhytI,p,q); % downsample from 490 hz to 512Hz

data = RhytI;
cfg.toi = linspace(10,48,38*EEG.srate);
temp = data;              
dat = [zeros(size(temp,1),10*fs),temp,zeros(size(temp,1),10*fs)]; 
padding = 2^nextpow2(length(dat));
cfg.pad = padding/fs;
options = {'pad', cfg.pad, 'padtype', cfg.padtype, 'freqoi', cfg.foi, 'polyorder', cfg.polyremoval};
[spectrum,foi,toi] = ft_specest_wavelet(dat, linspace(0,58,length(dat)), 'timeoi', cfg.toi, 'width', cfg.width, 'gwidth', cfg.gwidth,options{:});
spectrum(:,:,1:2*fs) = [];
spectrum(:,:,end-2*fs+1:end) = [];
RhytSpec{1} = spectrum;

[RhytII,fN] = audioread('\Quadruple Rhythm.wav');
RhytII      = RhytII(:,1);
RhytII      = envelope(RhytII,Wn,'analytic');
r1          = 10;
RhytII      = decimate(RhytII,r1,'fir'); % downsample from  44100 Hz to 4410 Hz
r2          = round(4410/512);
RhytII      = decimate(RhytII,r2,'fir'); % downsample from 4410Hz to 490Hz
p           = 256;
q           = 245;
RhytII      = resample(RhytII,p,q); % downsample from 490 hz to 512Hz

data = RhytII;
cfg.toi = linspace(10,46,36*EEG.srate);
temp = data;              
dat = [zeros(size(temp,1),10*fs),temp,zeros(size(temp,1),10*fs)]; 
padding = 2^nextpow2(length(dat));
cfg.pad = padding/fs;
options = {'pad', cfg.pad, 'padtype', cfg.padtype, 'freqoi', cfg.foi, 'polyorder', cfg.polyremoval};
[spectrum,foi,toi] = ft_specest_wavelet(dat, linspace(0,56,length(dat)), 'timeoi', cfg.toi, 'width', cfg.width, 'gwidth', cfg.gwidth,options{:});
spectrum(:,:,1:4*fs) = [];
spectrum(:,:,end-4*fs+1:end) = [];
RhytSpec{2} = spectrum;

%% Phase Calculation
freqbin1 = 5; % it changes based on the desired frequency
freqbin2 = 11; % it changes based on the desired frequency
for i=1:length(sub)
    Data{i,1} = squeeze(spec{i,1}(:,freqbin1,:));
    Data{i,2} = squeeze(spec{i,2}(:,freqbin2,:));
end
s1 = squeeze(RhytSpec{1}(:,freqbin1,:));
s1 = angle(s1);
s1 = s1';

s2 = squeeze(RhytSpec{2}(:,freqbin2,:));
s2 = angle(s2);
s2 = s2';
for i=1:length(sub)
    Data{i,1} = angle(Data{i,1});
    Data{i,2} = angle(Data{i,2});
end
for i=1:length(sub)
    for j=1:128
        pdata{i,1}(j,:) = circ_dist(Data{i,1}(j,:),s1);
        pdata{i,2}(j,:) = circ_dist(Data{i,2}(j,:),s2);
    end
end
for i=1:length(sub)
    ppdata{i,1} = exp(1i*pdata{i,1});
    ppdata{i,1} = mean(ppdata{i,1},2);
    ppdata{i,2} = exp(1i*pdata{i,2});
    ppdata{i,2} = mean(ppdata{i,2},2);
end

for i=1:length(sub)
    Phase{1}(:,i) = angle(ppdata{i,1});
    Phase{2}(:,i) = angle(ppdata{i,2});
end