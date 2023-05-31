function [Sniff,Sniff_smooth,frame_trigger_trial,frame_trigger,data,Sniff_time]=Read_Trial_Info(h5,path_used,pre,post,inh_detect,fps,wsFrameNumbers, pmtblank_dur)
%This function is to extract sniff trace of individual trials
%Modified for use in 2p imaging
%Sniff(trial x ms): sniff traces of each trial
% Sniff_smooth: smoothed version of sniff
% frame_trigger_trial : Frame trigger times for each trial
%frametrigger: frame triggers for tif stack (half the number of total
%frames in tiff stack because both green and red channel is recorded at
%each triggering)
%Sniff_time(trial x ms): time in data file for corresponding time bins in
%Sniff
%
% Mursel Karadas 2022

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

if ~exist('fps','var')
    % Fps is given or notS
    fps=30;
end
if ~exist('pmtblank','var')
    pmtblank_dur = 16; 
end
if ~exist('wsFrameNumbers','var')
    Blanked_recording = false;
else
    Blanked_recording = true;

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
fvOnTime=data.fvOnTime;
%%
%obtain frametrigger
record_onset=zeros(size(inh_onset));
frame_trigger= [];
frame_tol = 1e3/fps*1.2;  
delta_t = 1e3/fps;
frame_fv = zeros(num_trial,1);
for i=1:num_trial
    %i
    data_Events=h5read(h5,strcat(Keys{i},'/Events'));

    packet_sent_time=data_Events.packet_sent_time;
    sniff_samples=data_Events.sniff_samples;    %array for the duration of each packet
    record_onset(i)=packet_sent_time(1)-sniff_samples(1); %onset time for Sniff{trial(i)}
    trialframes = cell2mat(h5read(h5,strcat(Keys{i},'/frame_triggers')));
    if ~isempty(trialframes)
        frame_trigger = [frame_trigger; trialframes];
        [~,frame_fv(i)] = min(abs( frame_trigger -double(fvOnTime(i))));
    end
end
Blanked_Duration = 1e3;  % To identify each trial since voyeur start trial by FV opening.
%% Detecting missing frametrigger and fill the missing frametrigger using
%linear interploation
%some frametriggers are missed
%some other frametriggers are replaced by later time (usually 300 - 400ms
%later). This make the same frametrigger happing twice.
%Remove duplicate
if length(frame_trigger)~=length(unique(frame_trigger))
    fprintf('duplicated frametrigger \n')
end
frametrigger2 = sort(unique(frame_trigger));%If there is no duplicates, frametrigger2=frametrigger;
if Blanked_recording
    Frame_endindices = [find(diff(frametrigger2)>Blanked_Duration); length(frametrigger2)];
    if length(wsFrameNumbers) >1
        VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];
        fprintf(['VoyeurFrameTotal: ',num2str(sum(VoyeurFrameNumbers))]);
        fprintf(['WSFrameTotal: ',num2str(sum(wsFrameNumbers))]);
        [ProblematicTrials ,~,~]= find((VoyeurFrameNumbers(:)-wsFrameNumbers(:)) ~=0);
        if isempty(ProblematicTrials(:))
            fprintf('No frame loss in Voyeur \n');
        else
            fprintf(['Trials numbers :  ',num2str(ProblematicTrials'), '   have missing frames.','\n']);
            Frametrigger_temp = frametrigger2;
         for i = 1:length(ProblematicTrials)
          if ProblematicTrials(i)==1
                    ProblematicFrameBounds = [1 Frame_endindices(ProblematicTrials(i))];
          else
                ProblematicFrameBounds = [Frame_endindices(ProblematicTrials(i)-1)+1 Frame_endindices(ProblematicTrials(i))];
                frameidx=find(diff(Frametrigger_temp(ProblematicFrameBounds(1): ProblematicFrameBounds(2)))>frame_tol,1);
                 if ~isempty(frameidx)
                    frameidx =  frameidx + ProblematicFrameBounds(1) -1;
                 end
                fill_trial = 0;
                while ~isempty(frameidx)
                    fill=round((Frametrigger_temp(frameidx+1)-Frametrigger_temp(frameidx))/delta_t); %missing frames
                    errorvals = Frametrigger_temp(frameidx:frameidx+1);
                    Frametrigger_temp=[Frametrigger_temp(1:frameidx);...
                        round(Frametrigger_temp(frameidx)+delta_t:delta_t:Frametrigger_temp(frameidx)+delta_t*(fill-1))';...
                        Frametrigger_temp((frameidx+1):end)]; % make sure to add the right number of frames
                    fprintf(['We replaced ',num2str(errorvals'),' with ',...
                        num2str(Frametrigger_temp(frameidx:frameidx+fill)'),'\n']);
                    ProblematicFrameBounds  = [ProblematicFrameBounds(1) ProblematicFrameBounds(2)+fill-1];
                    frameidx=find(diff(Frametrigger_temp(ProblematicFrameBounds(1): ProblematicFrameBounds(2)))>frame_tol,1);
                    fill_trial = fill_trial + fill-1;
                    if ~isempty(frameidx)
                        frameidx =  frameidx + ProblematicFrameBounds(1) -1;
                    end
                    if length(errorvals') == length(Frametrigger_temp(frameidx:frameidx+fill)')
                        frameidx = [];
                    end
                end
                Frame_endindices = [find(diff(Frametrigger_temp)>1e3); length(Frametrigger_temp)];
                VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];
                if VoyeurFrameNumbers(ProblematicTrials(i)) == wsFrameNumbers(ProblematicTrials(i))
                   fprintf(['Trials numbers :  ',num2str(ProblematicTrials(i)'), '  matched to wavesurfer.','\n']);
                end
            end
            frametrigger2 = Frametrigger_temp;
    end
        end
        %% check again especially for first trial since Voyeur can drop many frames at the beginning
        if(VoyeurFrameNumbers(1)-wsFrameNumbers(1)) ~=0
            Nmiss = wsFrameNumbers(1) -VoyeurFrameNumbers(1);
            frametrigger2 = vertcat(linspace(frametrigger2(1)-delta_t-Nmiss*delta_t, frametrigger2(1)-delta_t, Nmiss)', frametrigger2);
            Frame_endindices = [find(diff(frametrigger2)>Blanked_Duration); length(frametrigger2)];
            VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];

        end

    else
        if ~sum( ~diff(frametrigger2)<1.2*1e3/fps &  ~diff(frametrigger2)>0.8*1e3/fps)
            fprintf('No problem in frames \n');
        end
    end
    
    %% Assign frame times to each trials by checking starttrial infos
    frame_trigger_trial = {};
    frames = frametrigger2;
    N_tiff = length(VoyeurFrameNumbers);
    data.Voyeur_frame_numbers = [];
    for i = 1: N_tiff
        Nf = VoyeurFrameNumbers(end-i+1);
        frames_i = frames(end-Nf+1:end);
        for j = 1:num_trial
            [val(j),frame_diff(j)] = min(abs(frames_i -double(data.starttrial(j))));
        end
        [~,trial_id] = min(val); 
        frame_trigger_trial{trial_id} = frames_i;
        frames(end-Nf+1:end) = [];
        data.Voyeur_frame_numbers(trial_id) = size(frames_i,1);
    end

else
    frame_trigger_trial = {};
    if nnz(diff(frametrigger2)>frame_tol)
        frameidx=find(diff(frametrigger2)>frame_tol,1);
        while ~isempty(frameidx)
            fill=round((frametrigger2(frameidx+1)-frametrigger2(frameidx))/delta_t); %missing frames
            errorvals = frametrigger2(frameidx:frameidx+1);
            frametrigger2=[frametrigger2(1:frameidx);...
                round(frametrigger2(frameidx)+delta_t:delta_t:frametrigger2(frameidx)+delta_t*(fill-1))';...
                frametrigger2((frameidx+1):end)]; % make sure to add the right number of frames
            fprintf(['We replaced ',num2str(errorvals'),' with ',...
                num2str(frametrigger2(frameidx:frameidx+fill)'),'\n']);
            %         frameidx=find(diff(frametrigger2)>35,1); Caused infinite loop in
            %         14393_190914_field4
            frameidx=find(diff(frametrigger2)>frame_tol,1);
        end
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

lost_packet = find(diff(packet_sent_time)~=sniff_samples(2:end))+1;
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
                for k=1:length(pos_lost)
                    j = pos_lost(k);
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
Sniff(Sniff_time==0)=NaN;
Sniff_time(Sniff_time==0) = NaN; %% Sniff_time =0 is problematic acqusition and cause wrong estiamtion for inh and fv bin, MK 21/07/10
row_NaN = isnan(mean(Sniff_time,2)); % We will add this to inh and fv bin;
locs_NaN = find(row_NaN ==1);
[inh_bin,~]=find(Sniff_time'==double(inh_onset)',length(inh_onset),'first');
inh_bin = insert_element(inh_bin,locs_NaN, zeros(length(locs_NaN),1));
% inh_bin=findfirst(Sniff_time==double(inh_onset),2);
[bin_fv,trial_id]=find(Sniff_time'==double(fvOnTime)',length(fvOnTime),'first');
% fv_bin=findfirst(Sniff_time==double(fvOnTime),2);
fv_bin = 1950*ones(length(fvOnTime),1);
fv_bin(trial_id) = bin_fv;
fv_bin = insert_element(fv_bin,locs_NaN, zeros(length(locs_NaN),1));

%% Find Stim trigger
stim_time = inh_onset + data.pulseOnsetDelay_1;
for i = 1:num_trial-1
    %stim_frame(i)=find(frame_trigger<stim_time(i+1), 1, 'last');%frame in which inh_onset is included
    %if data.amplitude_1(i) >0
        [data.stim_diff(i), ind] = min(abs(double(frame_trigger)+fps/2.0 -double(stim_time(i+1)- pmtblank_dur/2.0)));
        data.stim_frame(i) = ind;
    %end
end
%% Need to change this it causes problem in inh detection
inh_inSniff=inh_onset-record_onset(1);
Snifftemp=zeros(size(Sniff));
for i=1:num_trial
    if nnz(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1))>0
        
        Snifftemp(i,:)=sniff_all(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1));
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    end
end
inh_bin2= inh_bin;
%%
if inh_detect
    for kk = 2:num_trial
        kk
        respiratoryTrace = Snifftemp(kk,:);
        bmObj = breathmetrics(-1*respiratoryTrace', 1e3, 'rodentAirflow');
        bmObj.estimateAllFeatures(1,'simple', 0, 1);
        [~, index] = min(abs(bmObj.inhaleOnsets-2001));
        inh_bin2(kk) = bmObj.inhaleOnsets(index);
        %fig = bmObj.plotFeatures({'onset'});
    end

    data.inh_onset_voyeur=data.inh_onset;
    inh_onset=data.inh_onset+int32(inh_bin2(:)-inh_bin(:));
    data.inh_onset=inh_onset;
    
else
    data.inh_onset_voyeur=data.inh_onset;
end
%%
inh_inSniff=inh_onset-record_onset(1);
Sniff2=zeros(size(Sniff));
for i=2:num_trial
    if nnz(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1))>0
        
        Sniff2(i,:)=sniff_all(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1));
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    end
end
Sniff=single(Sniff2);
data.frametrigger = frame_trigger;
data.fv_bin = fv_bin;
Sniff_smooth = Sniff;
for i=2:num_trial
    bmObj = breathmetrics(-1*Sniff(i,:)', 1e3, 'rodentAirflow');
    bmObj.estimateAllFeatures(1,'simple', 0, 1);
    Sniff_smooth(i,:) = -1*bmObj.baselineCorrectedRespiration;
    [~, index] = min(abs(bmObj.inhaleOnsets-2001));
    data.pre_inhs{i} = bmObj.inhaleOnsets(1:index-1);
    data.post_inhs{i} = bmObj.inhaleOnsets(index:end);
end
cd(path_orig);


function [v] = insert_element(vin,locs, locsval)
   v = vin;
   for i =1: length(locs)
       v = vertcat(v(1:locs(i)-1), locsval(i), v(locs(i):end));
   end    
