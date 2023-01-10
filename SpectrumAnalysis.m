clear, clc

sub = 63:81;
sub([5,15]) = [];
FFT = cell(length(sub),2);
load('Chan128.mat'); % channel locations
chanlocs = cell(1,128);
for i=1:length(chanlocs)
    chanlocs{i} = chan(i).labels;
end 

%% Calculating the normalized FFT of each subject
for s=1:length(sub)
    address = ['\N',num2str(sub(s)),'\PreFFTdata.mat']; % Load preprocessed data of each subject
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
% Frequency Tagging approach  
    sigS001             = nanmean(sigS001,3);
    sigS001             = sigS001 - mean(sigS001,2);
    freqS001            = abs(fft(sigS001')/size(sigS001,2))';
    freqS001            = freqS001(:,1:size(sigS001,2)/2+1);
    freqS001(:,2:end-1) = 2*freqS001(:,2:end-1);
    fS001               = EEG.srate*(0:(size(sigS001,2)/2))/size(sigS001,2);       
    
    temp = freqS001;
    tmp  = find(fS001>=.5 & fS001<=40);
    tmp  = [tmp(1),tmp(end)];
    for i=tmp(1):tmp(2)
        freqS001(:,i) = temp(:,i) - (mean(temp(:,i-5:i-3),2) + mean(temp(:,i+3:i+5),2))/2;
    end
        
    sigS002             = nanmean(sigS002,3);
    sigS002             = sigS002 - mean(sigS002,2);
    freqS002            = abs(fft(sigS002')/size(sigS002,2))';
    freqS002            = freqS002(:,1:size(sigS002,2)/2+1);
    freqS002(:,2:end-1) = 2*freqS002(:,2:end-1);
    fS002               = EEG.srate*(0:(size(sigS002,2)/2))/size(sigS002,2);
    
    temp = freqS002;
    tmp  = find(fS002>=.5 & fS002<=40);
    tmp  = [tmp(1),tmp(end)];
    for i=tmp(1):tmp(2)
        freqS002(:,i) = temp(:,i) - (mean(temp(:,i-6:i-3),2) + mean(temp(:,i+3:i+6),2))/2;
    end
    
    FFT{s,1} = freqS001;
    FFT{s,2} = freqS002;
end

%% plot
cfg = [];
cfg.xlim = [0 6];
cfg.parameter = 'avg';
cfg.showlabels = 'yes';
cfg.showoutline = 'yes';

signalS001 = zeros(size(FFT{1,1},1),size(FFT{1,1},2),length(sub));
signalS002 = zeros(size(FFT{1,2},1),size(FFT{1,2},2),length(sub));

for s=1:length(sub)
    signalS001(:,:,s) = FFT{s,1};
    signalS002(:,:,s) = FFT{s,2};
end
temp1 = signalS001;
temp2 = signalS002;
signalS001 = nanmean(signalS001,3);
signalS002 = nanmean(signalS002,3);
    
EEG1 = pop_importdata('data',signalS001,'srate',512);
EEG1.chanlocs = chan;
fielddata1 = eeglab2fieldtrip(EEG1,'timelockanalysis','coord_transform');
fielddata1.time = fS001;
tmp = temp1;
num = zeros(size(tmp,1),1);
for q=1:size(tmp,1)
    num(q) = sum(~isnan(sum(tmp(q,:,:),2)),3);
end
fielddata1.std = nanstd(temp1(:,:,:),[],3)./repmat(sqrt(num),1,size(fielddata1.avg,2));
fielddata1.time = fS001;   

EEG2 = pop_importdata('data',signalS002,'srate',512);
EEG2.chanlocs = chan;
fielddata2 = eeglab2fieldtrip(EEG2,'timelockanalysis','coord_transform');
fielddata2.time = fS002;
tmp = temp2;
num = zeros(size(tmp,1),1);
for q=1:size(tmp,1)
    num(q) = sum(~isnan(sum(tmp(q,:,:),2)),3);
end
fielddata2.std = nanstd(temp2(:,:,:),[],3)./repmat(sqrt(num),1,size(fielddata2.avg,2));
fielddata2.time = fS002;   

temp = fielddata1.elec.pnt(:,1);
fielddata1.elec.pnt(:,1) = -fielddata1.elec.pnt(:,2);
fielddata1.elec.pnt(:,2) = temp;
fielddata2.elec.pnt = fielddata1.elec.pnt;
    
figure
ft_multiplotER(cfg,fielddata1);
figure
ft_multiplotER(cfg,fielddata2);