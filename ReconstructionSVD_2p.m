function [dff,trial_info,trials_read, meanAll,F]=ReconstructionSVD_2p(U,SV,roiMasks,trial_info,pre,post)
% This function reconstructs df/f using SVD compressed data. Mostly used
% for analyzing 2p imaging data
%
%Input
%U: spatial component obtained by SVD of raw imaging data
%SV: product of singular value and temporal component obtained by SVD of raw imaging data
%trial_info: information about imaging session
%roiMasks: binary mask of ROIs 2d(pixel x roi) or 3d matrix(pixel x pixel x roi)
%pre:time before inh_onset of each trial (ms)
%post:time after inh_onset of each trial (ms)
%
%Output
%dff: df/f (roi x time x trial)
%trial_info:structure containing trial information
%trials_read:0 or 1 to tell if the trials are correctly read
%meanAll: mean dff across all ROIs in each trial. Used for debugging
%
%Hirofumi Nakayama 2020

if ~exist('pre','var')
    %time before inh_onset of each trial (ms)
    pre=500;
end

if ~exist('post','var')
    %time after inh_onset of each trial (ms)
    post=1000;
end

inh_onset=trial_info.inh_onset;
frametrigger=trial_info.frametrigger;
num_trial=length(inh_onset);
%Temporal upsampling around inhalation onset, Sniff sampling rate is 1 kHz
up_sampling_fac=(1000/trial_info.fps); 

%roiMasks supposed to have shape of (pixel x roi)
%Reshape it if stack of binary images are given
if ndims(roiMasks)==3
   roiMasks = reshape(roiMasks,[],size(roiMasks,3)); 
end

%lower and upper singular values used
s_lower = 1;
s_upper =50;

Ucell=roiMasks'*U(:,s_lower:s_upper);

%number of frames pre- post- inhalation of each trial for upsampling
Nframe = size(SV,1);
pre_frame=floor(pre/trial_info.fps);
post_frame=floor(post/trial_info.fps);
dff=[];
meanAll=[];
F =zeros(size(Ucell,1),pre+post,num_trial);
for tr=1:num_trial
    inh=inh_onset(tr);
    inh_frame=find(frametrigger<inh, 1, 'last' );%frame in which inh_onset is included
    if ~isempty(inh_frame)&&inh_frame+post_frame-1<=Nframe&&inh_frame-pre_frame>=1
        sv_frame_range=inh_frame-pre_frame:inh_frame+post_frame-1;
        sv_upsample=interp1(SV(sv_frame_range,:),1:1/up_sampling_fac:(pre_frame+post_frame)); %upsample imaging data
        
        fcell=sv_upsample(single([inh-pre:inh+post-1])-frametrigger(inh_frame-pre_frame),(s_lower:s_upper))*Ucell';
        baseline_frame=1:pre; %
        %dff(cellID, time from inh_onset(ms), trial)
        dff(:,:,tr)=((fcell-mean(fcell(baseline_frame,:)))./mean(fcell(baseline_frame,:)))';
        F(:,:,tr) = fcell';
        trials_read(tr)=true;
        meanAll(:,tr) = (mean(fcell(baseline_frame,:)));
    else
        dff(:,:,tr)=zeros(size(Ucell,1),pre+post);
        F(:,:,tr)=zeros(size(Ucell,1),pre+post);
        fprintf('trial %d inh_frame was not included in tiff stack',tr)
        trials_read(tr)=false;
        meanAll(:,tr) = zeros(size(Ucell,1),1);
    end
    
end


