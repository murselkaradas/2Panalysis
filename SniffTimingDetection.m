function FS=SniffTimingDetection(Sniff,varargin)
% Extract timing of events in sniff signal
% Input
%     Sniff: (trial x time) time include -1000 to +2000ms from inhalation onset
% Output
%     FS: Structure containing timings of different events
%     
% Hirofumi Nakayama 2019

if numel(varargin)==0
    trial=1:size(Sniff,1);
    fig=0;
elseif numel(varargin)==1
    trial=varargin{1};
    fig=0;
elseif numel(varargin)==2
    [trial,fig]=varargin{:};
end

%This part needs to be fixed
s_var=var(Sniff);
%Detect the timing of inhalation trigger
LB = 300; %ms, lower bound for inhalation trigger search
UB = 2000;%ms, Upper bound for inhalation trigger search
[~,inh_trig]=min(s_var(LB+1:UB));
inh_trig=inh_trig+LB;
pre=inh_trig-1;
post=size(Sniff,2)-inh_trig;

trials_unread=max(Sniff,[],2)==0;
trial(ismember(trial,find(trials_unread)))=[];

SF=zeros(size(Sniff));
for i=trial
    %detrend first to prevent linear y-shift
    %Modified 10/11/2019
    ds = detrend(Sniff(i,:));
    trend = Sniff(i,:)-ds;
    %Some rigs have discontinuous trend due to needle valves
    %Use median filtered signal to subtract discontinous baseline.
    tmp = butterworthed_signal(ds,1000,2,[0.5,30])+ trend;
    %Corrected 012520 medfilt1(tmp - medfilt1(tmp,750),30) loose sharp peaks
    % SF(i,:)=medfilt1(tmp - medfilt1(tmp,750),30);%Medfilt40 to remove random peaks   
    SF(i,:)=medfilt1(tmp - medfilt1(tmp,750),5);
end

SF_d1=SF(:,2:end)-SF(:,1:end-1);SF_d1=[SF_d1(:,1),SF_d1];
SF_d2=SF_d1(:,2:end)-SF_d1(:,1:end-1);SF_d2=[SF_d2(:,1),SF_d2];

%Todo: set starting point based on histogram
th_start=-500;
Thres = cell(size(Sniff,1),1);
Position = cell(size(Sniff,1),1);
Inh_onset = cell(size(Sniff,1),1);
Inh_offset = cell(size(Sniff,1),1);
Inh_dur = cell(size(Sniff,1),1);
Sniff_dur = cell(size(Sniff,1),1);
Peaks_inh = cell(size(Sniff,1),1);
Peaks_exh = cell(size(Sniff,1),1);
Ind_inhonset = zeros(size(Sniff,1),1);

FS.inh_onset=zeros(size(Sniff,1),1);

Peak_threshold = 0.5; % normalized sniff threshold for peak detection
Min_Peak_Distance = 80; %ms, corresponds to maximum sniff freq.
for tr=trial
    clear position inh_onset inh_offset inh_dur sniff_dur scale_d
    
    %worked in 1p
    %     [~,~,peaks_exh,peaks_inh] = peakdet(SF(tr,:)-mean(SF(tr,:)), 100, 'th', 40);
    %worked in 2p
    %sf = (SF(tr,:)-mean(SF(tr,:)))/std(SF(tr,:));
    sf = normalize(SF(tr,:)); % normalizing the data to have mean 0 and standard deviation 1
    [~,~,peaks_exh,peaks_inh] = findpeak_pos_neg(sf, Peak_threshold,Min_Peak_Distance);
    if ~isempty(peaks_exh)&&~isempty(peaks_inh)
        if min(peaks_exh)>990
            peaks_exh=[960;peaks_exh];
        end
%         
        %Keep this part of code in case below modification screw up sniff
        %detection in 1p rig
        %         [~,~,peaks_pos_d1,peaks_neg_d1] = peakdet(SF_d1(tr,:)-mean(SF_d1(tr,:)), 5, 'th', 40); %threshold value 5 may need to be adjusted
        
        %Normalize SF_d1 by the amplitude of largest negative peaks
        %set threshold to 0.5
%         sfd = SF_d1(tr,:)-mean(SF_d1(tr,:));
%         for i=1:length(SF_d1(tr,:))-1000
%             scale_d(i+500) = std(SF_d1(tr,i:i+1000));
%         end
%         scale_d(1:500) = scale_d(501);
%         scale_d(end+1:end+500) = scale_d(end);
%         sfd = sfd./scale_d;
%         scale = std(SF(tr,:));
        
        sfd1 = normalize(SF_d1(tr,:));
        %abs(SF(tr,:)/scale) to detect peaks around SF=0
        %medfilt1(SF_d1(tr,:)./scale_d,200) to avoid random peaks in
        %pleatue
        
        %[~,~,peaks_pos_d1,peaks_neg_d1] = peakdet(SF_d1(tr,:)./scale_d + abs(SF(tr,:)/scale)+medfilt1(SF_d1(tr,:)./scale_d,200), 1, 'th', 40);%Threshold value is fraction to std
        [~,~,peaks_pos_d1,peaks_neg_d1] = findpeak_pos_neg(sfd1, 1,Min_Peak_Distance);

        %So far, use zero-crossing as the inhalation onset/offset
        
        %SF(tr,peaks_onset)
        %set adaptive threshold so that area below and above threshold becomes same
        
        %                     figure;plot(SF(tr,:));
        %                 hold on;plot(Sniff_filt_d1(tr,:)*10)
        
        
        %inh_peaks are mostly accurate
        %sometimes there are multiple exhalation peaks in a sniff
        %in some sniff, inh_peak are lost,
        %sniff1: inh, exh
        %sniff2: exh
        %simple zero crossing can't detect above errors
        
        %Some times FV trigger is detected by peaks_inh (this occur before peaks_exh)
        %skip such inhaltion
        %peaks_inh(1)and inh(end) are ignored because it may be very close to the edge and
        %may not have inh_onset/inh_offset
        %of the sniff.
        th=mean(SF(tr,:));
        
        exh1=zeros(length(peaks_inh),1);
        ind_exclude=false(length(peaks_inh),1);
        %ind_looked=find(peaks_inh(1:end-1)<2500)';
        ind_looked=find(peaks_inh(1:end)<size(Sniff,2))';
        ind_looked=ind_looked(2:end);
        %for i=2:length(peaks_inh)-1
        %Calculate exh1 exhs_1sniff peaks_inh ind_exclude
        for i=2:length(ind_looked)
            if ~isempty(find(peaks_exh<peaks_inh(i)))
                exh1(i)= peaks_exh(find(peaks_exh<peaks_inh(i),1,'last'));
            else
                exh1(i)=0;
            end
            if i==max(ind_looked)
                exhs_1sniff=peaks_exh((peaks_exh>peaks_inh(i)&peaks_exh<size(SF,2)))-exh1(i);
                
            else
                exhs_1sniff=peaks_exh((peaks_exh>peaks_inh(i)&peaks_exh<peaks_inh(i+1)))-exh1(i);
                
            end
            
            if peaks_inh(i)<1500&&i<max(ind_looked)%remove peaks_inh that is supposed to be FV onset
                
                if isempty(exhs_1sniff)&&(abs(SF(tr,peaks_inh(i+1)))<0.25*abs(mean(abs(SF(tr,peaks_inh)))))
                    ind_exclude(i+1)=true;%remove peaks_inh that is supposed to be FV onset
                end
            end
        end
        
        %In case mice hold a sniff and no inh_peak in 0-1000ms look for
        %exh_peak before 1000ms
        try                                                                         
            if peaks_inh(1)>1000
                i=1;
                if ~isempty(find(peaks_exh<peaks_inh(i)))
                    exh1(i)= peaks_exh(find(peaks_exh<peaks_inh(i),1,'last'));
                else
                    exh1(i)=0;
                end
                exhs_1sniff=peaks_exh((peaks_exh>peaks_inh(i)&peaks_exh<peaks_inh(i+1)))-exh1(i);
                if isempty(exhs_1sniff)&&(abs(SF(tr,peaks_inh(i+1)))<0.25*abs(mean(abs(SF(tr,peaks_inh)))))
                    ind_exclude(i+1)=true;%remove peaks_inh that is supposed to be FV onset
                end
            end
        catch
            %        figure;plot(SF(tr,:))
            %        a=1;
        end
        exh1(ind_exclude)=[];
        peaks_inh(ind_exclude)=[];
        
        %Check the duplicate of exh1. This indicate FV trigger is still
        %classified as inhalation
        ind_duplicate=(exh1-[100000;exh1(1:end-1)])==0;
        exh1(ind_duplicate)=[];
        peaks_inh(ind_duplicate)=[];
        
        if length(exh1)>1
            indices=find(exh1(1:end-1)>0)';
            %         else
            %             indices=1:1;
            %         end
            
            
            for i=indices
                th_val=200;
                clear pos
                clear area
                
                %use fixed baseline 0 as threshold
                thres(i)=th;
                try
                    if length(peaks_inh)==1
                        sn=SF(tr,peaks_inh(i):end);
                        peaks_inh(2)=size(SF,2);
                    else
                        sn=SF(tr,peaks_inh(i):peaks_inh(i+1));
                    end
                    exhs_1sniff=peaks_exh((peaks_exh>peaks_inh(i)&peaks_exh<peaks_inh(i+1)))-exh1(i);
                catch
                    %                 ffprintf('error detecting peaks_inh in trial %d',d)
                end
                sn2=SF(tr,exh1(i):peaks_inh(i+1))-thres(i);
                pn=sn2>0;
                
                tcross=find(pn(1:end-1)-pn(2:end));
                
                
                if exh1(i)+tcross(1)>peaks_inh(i)%first threshold-cross is missing
                    tcross=[0,tcross];
                end
                
                if length(tcross)==1
                    try
                        sn2=SF(tr,exh1(i):end)-thres(i);
                        pn=sn2>0;
                        tcross=find(pn(1:end-1)-pn(2:end));
                        tcross=tcross(1:3);
                    catch
                        
                        tcross(2)=length(pn)-1;
                        tcross(3)=length(pn);
                    end
                    
                end
                
                %In case multiple exhalation peaks detected in a single sniff
                if length(exhs_1sniff)>=2
                    
                    %Distance between multiple sniff
                    edist=diff(exhs_1sniff);
                    
                    %First exh_peaks before >100 interval
                    %Corresponding to multiple exh peaks in a exhalation cycle
                    exh_multi_ind=find(edist>100,1);
                    %if ~isempty(exh_multi_ind)
                    %skip this part if threshold crossing is correctly determined
                    if ~isempty(exh_multi_ind)&&length(tcross)>=4
                        %pnd:negative peak of 1st derivative after
                        %exhs_1sniff(exh_multi_ind)+exh1
                        try
                            neg1d_after_exh=peaks_neg_d1(peaks_neg_d1>exhs_1sniff(exh_multi_ind)+exh1(i));
                            %near_neg_peak_1d=pnd(find(pnd>exhs_1sniff(exh_multi_ind),1));
                            tcross(3)=find(SF_d1(tr,neg1d_after_exh(1):end)>0,1)+neg1d_after_exh(1)-exh1(i);
                        catch
                            %                        ffprintf('Trial%d sniff%d exh1(i)=%d: error in tcross(3) = xxx\n',tr,i,exh1(i))
                        end
                    end
                    
                elseif isempty(exhs_1sniff)
                    %No exhalation in a inter-inhalation interval
                    %due to FV trigger detected as inhalation
                    %in case exh is too weak, threshold crossing should still happen
                end
                
                try
                    
                    position(i,:)=exh1(i)+[tcross(1)+1,tcross(2),tcross(3)];
                    inh_onset(i)=exh1(i)+tcross(1)+1;
                    inh_offset(i)=exh1(i)+tcross(2);
                    inh_dur(i)=tcross(2)-tcross(1);
                    sniff_dur(i)=tcross(3)-tcross(1);
                catch
                    %                 ffprintf('error detecting inh_exh in trial%d',tr)
                end
                %Define baseline from first and second order derivative
                %first local mimimun of the first derivative after inh_peaks
                %(or second zero crossing of the second derivatives after inh_peaks)
            end
        end
        %--------------------------------------------------------
        %detect index for the sniff for inh_trigger
        %In case onset transient of FV is detected as inh
        %Chose the next sniff
        %[~,ind]=min(abs(inh_onset-inh_trig));
        
        ind_first=find((peaks_inh-inh_trig)>0,1);
        %[~,ind_first]=min(abs(peaks_inh-inh_trig));%Why this was used?
        %--------------------------------------------------------
        
        try
            Ind_inhonset(tr)=ind_first;
            Peaks_inh{tr}=peaks_inh;
            Peaks_exh{tr}=peaks_exh;
            Thres{tr}=thres;
            Position{tr}=position;
            Inh_onset{tr}=inh_onset;
            Inh_offset{tr}=inh_offset;
            Inh_dur{tr}=inh_dur;
            Sniff_dur{tr}=sniff_dur;
        catch
            %             ind_first
        end
        
        try
            
            FS.position(tr,:)=position(ind_first,:)';
            FS.peak_inh(tr)=peaks_inh(ind_first);
            FS.peak_inh2(tr)=peaks_inh(ind_first+1);
            FS.peak_exh(tr)=peaks_exh(ind_first);
            FS.peak_exh2(tr)=peaks_exh(ind_first+1);
            FS.thres(tr)=thres(ind_first);
            FS.inh_onset(tr)=inh_onset(ind_first);
            FS.inh_onset2(tr)=inh_onset(ind_first+1);
            FS.inh_offset(tr)=inh_offset(ind_first);
            FS.inh_offset2(tr)=inh_offset(ind_first+1);
            FS.inh_dur(tr)=inh_dur(ind_first);
            FS.inh_dur2(tr)=inh_dur(ind_first+1);
            FS.sniff_dur(tr)=sniff_dur(ind_first);
            FS.sniff_dur2(tr)=sniff_dur(ind_first+1);
            FS.inh_trig(tr)=inh_trig;
            FS.first_sniff{tr}=SF(tr,position(ind_first,1)-10:position(ind_first,3));
            FS.sniff{tr} = SF(tr,:);
        catch
            fprintf('error in calculating sniff features trial %d',tr)
        end
        
    else
        
    end
end
FS.Inh_onsets_all=Inh_onset;
FS.Inh_offsets_all=Inh_offset;
FS.Inh_peaks_all = Peaks_inh;
FS.Exh_peaks_all = Peaks_exh;



end
function filtsig = butterworthed_signal(inputsig,fs,order,band)
    [b,a]=butter(order,band./(fs/2),'bandpass');
    filtsig=filtfilt(b,a,double(inputsig));  %filtered signal
end

function [pos_pks,neg_pks,pos_id,neg_id] = findpeak_pos_neg(sig, t, mdist)
    [pos_pks,pos_id] = findpeaks(sig,'MinPeakHeight',t,'MinPeakDistance',mdist);
    [neg_pks,neg_id] = findpeaks(-sig,'MinPeakHeight',t,'MinPeakDistance',mdist);
end