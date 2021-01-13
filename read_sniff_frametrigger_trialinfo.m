function [Sniff,frame_trigger,data,Sniff_time]=read_sniff_frametrigger_trialinfo(h5,path_used,pre,post,inh_detect)
%This function is to extract sniff trace of individual trials
%Modified for use in 2p imaging
%Sniff(trial x ms): sniff traces of each trial
%frametrigger: frame triggers for tif stack (half the number of total
%frames in tiff stack because both green and red channel is recorded at
%each triggering)
%Sniff_time(trial x ms): time in data file for corresponding time bins in
%Sniff
%
%Hirofumi Nakayama 2020

%Duration before inahlation to be analyzed
if ~exist('pre','var')
    pre=1000;
end

%Duration after inhalation to be analyzed
if ~exist('post','var')
    post=1000;
end

if ~exist('inh_detect','var')
    %Whether or not detect inhalation by threshold crossing
    inh_detect=true;
end

path_orig=pwd;
cd(path_used)

try
    data=h5read(h5,'/Trials');
catch
    error('%s not found in %s',h5, path_used)
end
%Need to check if this manipulation is ok for 2p data
inh_onset=double(data.inh_onset);
num_trial=length(inh_onset);

h_info=h5info(h5);
h_info=h_info.Groups;
% Keys is Number_of_trials x 1 
Keys={h_info.Name};%h5 file from 2p rig doesn't record trial0 analog signal
%%
%obtain frametrigger
record_onset=zeros(size(inh_onset));
frame_trigger= [];
for i=1:num_trial
    %i
    data_Events=h5read(h5,strcat(Keys{i},'/Events'));
    
    packet_sent_time=data_Events.packet_sent_time;
    sniff_samples=data_Events.sniff_samples;    %array for the duration of each packet
    record_onset(i)=packet_sent_time(1)-sniff_samples(1); %onset time for Sniff{trial(i)}
    frame_trigger = [frame_trigger; cell2mat(h5read(h5,strcat(Keys{i},'/frame_triggers')))];
end

%% Detecting missing frametrigger and fill the missing frametrigger using
%linear interploation
%some frametriggers are missed
%some other frametriggers are replaced by later time (usually 300 - 400ms
%later). This make the same frametrigger happing twice.

%Remove duplicate
if length(frame_trigger)~=length(unique(frame_trigger))
    fprintf('duplicated frametrigger')
end
frametrigger2 = sort(unique(frame_trigger));%If there is no duplicates, frametrigger2=frametrigger;

% frame_tol = 34
frame_tol = 40;  % what is 40 
if nnz(diff(frametrigger2)>frame_tol)
    frameidx=find(diff(frametrigger2)>frame_tol,1);
    while ~isempty(frameidx)
        fill=round((frametrigger2(frameidx+1)-frametrigger2(frameidx))/33.3);
        errorvals = frametrigger2(frameidx:frameidx+1);
        frametrigger2=[frametrigger2(1:frameidx);...
            round(frametrigger2(frameidx)+33.3:33.3:frametrigger2(frameidx)+33.4*(fill-1))';...
            frametrigger2((frameidx+1):end)]; % make sure to add the right number of frames
        fprintf(['We replaced ',num2str(errorvals'),' with ',...
            num2str(frametrigger2(frameidx:frameidx+fill)'),'\n']);
        %         frameidx=find(diff(frametrigger2)>35,1); Caused infinite loop in
        %         14393_190914_field4
        frameidx=find(diff(frametrigger2)>frame_tol,1);
    end
end

frame_trigger = frametrigger2;
%%
%Inhalation timing within sniff traces of each trial
tot_pre_post = pre + post;
inh_onset_local=inh_onset-record_onset;
Sniff=zeros(num_trial,pre+post);
Sniff_time=zeros(num_trial,pre+post);

%check missing sniff packets
packet_sent_time=[];
sniff_samples=[];
trial_index = [];
trial_subind = [];
for i=1:num_trial
    %check missing sniff packet
    events = h5read(h5,strcat(Keys{i},'/Events'));
    packet_sent_time = [packet_sent_time;events.packet_sent_time];
    sniff_samples = [sniff_samples;events.sniff_samples];
    trial_index = [trial_index;i*ones(length(events.sniff_samples),1)];
    trial_subind = [trial_subind;[1:length(events.sniff_samples)]'];
end

lost_packet = find(diff(packet_sent_time)~=sniff_samples(2:end))+1; %why  +1 
sniff_all=[];
for i=1:num_trial
    st=zeros(1,pre+post);
    try
        sniffcell=h5read(h5,strcat(Keys{i},'/sniff'));
        if ~isempty(lost_packet)
            pos_lost=lost_packet(trial_index(lost_packet)==i);%position of lost packet in time from start
            
            %             subind=trial_subind(lost_packet);
            %             subind(trial_index(lost_packet)~=i)=[];
            if ~isempty(pos_lost)
                for j=pos_lost
                    fprintf('trial %d, %d-th packet lost',i,trial_subind(j))
                    sniffcell{trial_subind(j)}=[int16(zeros((packet_sent_time(j)-packet_sent_time(j-1))-length(sniffcell{trial_subind(j)}),1));sniffcell{trial_subind(j)}];
                end
                
            end
        end
        sniff0=cell2mat(sniffcell);
        sniff_all=[sniff_all;sniff0];
        %pad with 10000 zeros before and after sniff
        sniff=[zeros(2*tot_pre_post,1);sniff0;zeros(2*tot_pre_post,1)];
        sniff_pos=[(-2*tot_pre_post+1):0,1:length(sniff),length(sniff)+1:length(sniff)+(2*tot_pre_post)];
        
        sniff_range=inh_onset_local(i)-pre:inh_onset_local(i)+post-1;
        
        st=sniff(ismember(sniff_pos,sniff_range));
        
        Sniff(i,:)=st;
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    catch
        sprintf('Error in trial %d',i)
    end
    
end

% what is -29999  and 30000
sniffall_time=[(-5*tot_pre_post+1):0,1:length(sniff_all),length(sniff_all)+1:length(sniff_all)+(5*tot_pre_post)];
sniff_all=[int16(zeros((5*tot_pre_post),1));sniff_all;int16(zeros((5*tot_pre_post),1))];

%%
%find real inhalation onset based on the threshold crossing of sniff signal
num_trial=size(Sniff,1);
inh_onset=data.inh_onset;
fvOnTime=data.fvOnTime;

%Check sniff alignment to fvOnTime or inh_onset
inh_bin=findfirst(Sniff_time==double(inh_onset),2);
fv_bin=findfirst(Sniff_time==double(fvOnTime),2);

offset=mean(Sniff(:));
pn=Sniff>offset;
Nfilt =30;
filt=[ones(1,Nfilt),-1*ones(1,Nfilt)]; % why 30, small cause false detection,  hand tuned 30 ms
trig=zeros(size(Sniff));
for t=1:size(Sniff,2)-2*Nfilt
    for tr=1:size(Sniff,1)
        trig(tr,t+Nfilt)=pn(tr,t:t+2*Nfilt-1)*filt';
        trig(tr,1:fv_bin(tr)+10)=0;
    end
end
inh_bin2=findfirst(trig==30,2);

if inh_detect
    data.inh_onset_voyeur=data.inh_onset;
    inh_onset=data.inh_onset+int32(inh_bin2-inh_bin);
    data.inh_onset=inh_onset;
    
    inh_inSniff=inh_onset-record_onset(1);
    Sniff2=zeros(size(Sniff));
    for i=2:num_trial
        if nnz(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1))>0
            Sniff2(i,:)=sniff_all(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1));
            Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
            
        end
    end
    Sniff=single(Sniff2);
end


data.frametrigger = frame_trigger;
cd(path_orig);
