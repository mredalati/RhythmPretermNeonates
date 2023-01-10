%% Load the data
load('data.mat'); % Load the data for each subject (channels*samples)
load('Chan128.mat') % Load the channel locations and labels
EEG = pop_importdata('data',data,'setname','Signal','srate',1000);
EEG.chanlocs = chan;
EEG.event = event;
[EEG.event.latency] = EEG.event.sample;
EEG.event = rmfield(EEG.event,'sample');
% There are two types of events (Duple/Triple Rhythm (S001) and Quadruple
% Rhythm (S002). S003 is trigger for silent periods between blocks.
for i=length(EEG.event):-1:1   
    if (strcmp(EEG.event(i).type,'S001') || strcmp(EEG.event(i).type,'S002') || ...
            strcmp(EEG.event(i).type,'S003'))
    else
        EEG.data(:,EEG.event(i+1).latency-1:-1:EEG.event(i).latency) = [];
        temp = length(EEG.event(i+1).latency-1:-1:EEG.event(i).latency);
        for j=i+1:length(EEG.event)        
            EEG.event(j).latency = EEG.event(j).latency - temp;                    
        end
        EEG.pnts = size(EEG.data,2);
        EEG.times = linspace(0,round(EEG.pnts/1000),EEG.pnts);
        EEG.event(i) = [];
    end    
end

%%
chanlocs = cell(1,EEG.nbchan);
for i=1:length(chanlocs)
    chanlocs{i} = EEG.chanlocs(i).labels;
end

channel  = [14,17,21,43,48,38,44,49,56,107,113,120,1,32,25,8,121,114,...
    99,94,88,81,73,68,63,119,125:128]; % Outer ring channels
E        = 'E';
chan     = cell(length(channel),1);
for i=1:length(chan)
    chan{i} = [E,num2str(channel(i))];
end
[~,loc]  = ismember(chan,chanlocs);
loc(loc==0) = [];
EEG.data(loc,:) = [];
EEG.chanlocs(loc) = [];
EEG.nbchan = EEG.nbchan - length(loc);

EEG = pop_eegfiltnew(EEG,.5,0);
EEG = pop_eegfiltnew(EEG,0,45);
EEG = pop_eegfiltnew(EEG,48,52,[],1);
EEG = pop_resample(EEG,512);

for i=1:length(EEG.event)
    EEG.event(i).latency = round(EEG.event(i).latency);
end

%% Artifact rejection
theta = 100;

trialmean = zeros(EEG.nbchan,length(EEG.event));
for i=1:length(EEG.event)
    if i~=length(EEG.event)
        trialmean(:,i) = mean(abs(EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency-1)),2);
    else
        trialmean(:,i) = mean(abs(EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+floor(38*EEG.srate))),2);
    end
end
rejTrial = trialmean>30;
for i=1:length(EEG.event)
    if sum(rejTrial(:,i))>=EEG.nbchan/2
        rejTrial(:,i) = 1;
    else
        rejTrial(:,i) = 0;
    end
end
rejTrial = mean(rejTrial);

% Applying AB algorithm and remove artifacts
for i=1:length(EEG.event)
    if i==length(EEG.event) && strcmp(EEG.event(end).type,'S001') % Duple/Triple Rhythm
        if rejTrial(i)
            EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(38*EEG.srate)) = NaN;
        else
            x = EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(38*EEG.srate));
            EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(38*EEG.srate)) = AB(x,theta);                     
        end
    elseif i==length(EEG.event) && strcmp(EEG.event(end).type,'S002') % Quadruple Rhythm
        if rejTrial(i)
            EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(36*EEG.srate)) = NaN;
        else
            x = EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(36*EEG.srate));
            EEG.data(:,EEG.event(i).latency:EEG.event(i).latency+...
                floor(36*EEG.srate)) = AB(x,theta);    
        end
    else
        if rejTrial(i)
            EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency-1) = NaN;
        else
            x = EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency-1);
            EEG.data(:,EEG.event(i).latency:EEG.event(i+1).latency-1) = AB(x,theta);
        end
    end    
end

% Average Reference
channel = 1:EEG.nbchan;
sig = EEG.data;
for i=1:EEG.nbchan
    temp = sig(channel~=i,:);
    temp = mean(temp);
    EEG.data(i,:) = sig(i,:) - temp;        
end