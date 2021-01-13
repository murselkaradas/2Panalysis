clear all

%% Add backend functions 
addpath(genpath('C:\Users\murse\Documents\MATLAB\2P_Analysis'))  

sn='session_analysis_062119.m';
%{
Odors used
Ethylbutyrate 1e-8M (4.86ul in 5.5ml for 1e-7M solution)
PropionicAcid 1e-8M (17.6ul in 5.5ml for 1e-7M solution)
ButyricAcid 4.5e-10M (5.5ul in 5.5ml for 4.5e-9M solution)
Benzaldehyde 4.9e-10M (2.2ul in 5.5ml for 4.9e-9M solution)
Methylvalerate 1e-9M (6.7ul in 5.5ml for 1e-7M solution)
Ethyltiglate 1e-9M (18.4ul in 5.5ml for 1e-7M solution)

Use longer ITI 12sec
Morphing stimuli
Ethylbutyrate 1e-8M (4.86ul in 5.5ml) and PropionicAcid 1e-8M (17.6ul in 5.5ml)
ButyricAcid 4.5e-10M (5.5ul in 5.5ml)and Benzaldehyde 4.9e-10M (2.2ul in 5.5ml)
Empty 100,90,71,50,29,10,0,0,10,29,50,71,90,100

PropionicAcid -> Ethylbutyrate
Benzaldehyde -> ButyricAcid
Ethylbutyrate -> PropionicAcid
ButyricAcid -> Benzaldehyde
Repeat this block twice / session
%}



direct=configPath();
%%
%Reconstruct all signals using compressed data
%Data compsression done on cluster

path_h5='C:\Users\hnaka\Dropbox\MATLAB\2P\2P_data\HN1953\190617';
%%
%Parameters to be changed
h5_name_PID='1953_01_D2019_6_17T11_18_1_PID_Morphing.h5';
field=1;

switch field
    case 1
        %field 1 morphing EB/PPA, Benzaldehyde/Butyric acid
        layer='mt';
        load('HN1953_190617_field1_svd.mat','U','SV','svals')
        clear H5
        H5{1}='1953_1_01_D2019_6_17T12_27_16_odor.h5';
        H5{2}='1953_1_03_D2019_6_17T13_9_7_odor.h5';
        sess_used=1:numel(H5);
        % h5_name='1953_1_01_D2019_5_17T11_38_22_odor.h5';
        name_data='1953_190617_field1_EB-PPA_ButAcd-BzAld_';
        % [cellMask,cellMask_sep]=create_ROImask_manual('HN1953_190603_field1_00003_ROIset');%Remove .tif from name
        [cellMask,cellMask_sep]=create_ROImask_manual('HN1953_190617_field1_RoiSet');%Remove .tif from name
        %}
        
    case 2
        %field2 glomerular layer of fiel1
        layer='glom';
        load('HN1953_190617_field2_svd.mat','U','SV','svals')
        clear H5
        H5{1}='1953_1_02_D2019_6_17T12_49_42_odor.h5';
        H5{2}='1953_1_04_D2019_6_17T13_27_16_odor.h5';
        H5{3}='1953_1_01_D2019_6_17T13_43_15_odor.h5';
        sess_used=1:2;
        % h5_name='1953_1_01_D2019_5_17T11_38_22_odor.h5';
        name_data='1953_190617_field2_EB_PPA_ButAcd-BzAld_';
        [cellMask,cellMask_sep]=create_ROImask_manual('HN1953_190617_field2_RoiSet');%Remove .tif from name
        
        
    case 3
        %field2 - MVT / ET habituation session
        layer='glom';
        load('HN1953_190617_field2_svd.mat','U','SV','svals')
        clear H5
        H5{1}='1953_1_02_D2019_6_17T12_49_42_odor.h5';
        H5{2}='1953_1_04_D2019_6_17T13_27_16_odor.h5';
        H5{3}='1953_1_01_D2019_6_17T13_43_15_odor.h5';
        sess_used=3:3;
        name_data='1953_190617_field2_MVT_ET_habituation_';
        [cellMask,cellMask_sep]=create_ROImask_manual('HN1953_190617_field2_RoiSet');%Remove .tif from name
otherwise
        disp('field not exists')
end

% [M,m,frameidx]=checkMissingData_2p()
%%


cd(path_h5)
Names = dir([strcat('*.h5')]);
Names={Names.name}';

% odor_set={'empty','EB100_PPA0','EB90_PPA10','EB71_PPA29','EB50_PPA50','EB29_PPA71','EB10_PPA90','EB0_PPA100'};
odor_set={'empty','EB0_PPA100','EB10_PPA90','EB29_PPA71','EB50_PPA50','EB71_PPA29','EB90_PPA10','EB100_PPA0','But0_Bz100','But10_Bz90','But29_Bz71','But50_Bz50','But71_Bz29','But90_Bz10','But100_Bz0'};
odor_num=[1,repmat([[1:8],[1,9:15],[1,fliplr(2:8)],[1,fliplr(9:15)]],1,2),1];

%cellMask_vec: n-th colum for cell n, 1/a in pixels, 0 otherwise
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);

pre=500;
post=1500;
Sniff={};

% for sess=sess_used
%     %     h5_name=Names{contains(Names,sprintf('1953_1_0%d',sess))};
%     h5_name=H5{sess};
%     %Calculate df/f for the first 1sec, 2sec, ...
%     [data,path_used]=read_h5(h5_name,path_h5);
%     
%     [sniff{sess,1},frametrigger,sniff_time{sess,1}]=Sniff_frametrigger_Extraction(h5_name,path_used,pre,post);
%     
%     num_trial=size(sniff{sess},1);
%     inh_onset=data.inh_onset;
%     fvOnTime=data.fvOnTime;
%     % od=data.olfas0x3Aolfa_00x3Aodor';
%     
%     odorconc=data.odorconc;
%     
%     mean_inh_onset=mean(findfirst(sniff_time{sess}==double(inh_onset),2));
%     mean_fvOnTime=mean(findfirst(sniff_time{sess}==double(fvOnTime),2));
%     
%     %Check sniff alignment to fvOnTime or inh_onset
%     inh_bin=findfirst(sniff_time{sess}==double(inh_onset),2);
%     fv_bin=findfirst(sniff_time{sess}==double(fvOnTime),2);
%     
%     offset=mean(sniff{sess}(:));
%     pn=sniff{sess}>offset;
%     filt=[ones(1,30),-1*ones(1,30)];
%     trig=zeros(size(sniff{sess}));
%     for t=1:size(sniff{sess},2)-length(filt)
%         for tr=1:size(sniff{sess},1)
%             trig(tr,t+length(filt)/2)=pn(tr,t:t+length(filt)-1)*filt';
%             trig(tr,1:fv_bin(tr)+10)=0;
%         end
%     end
%     inh_bin2=findfirst(trig==30,2);
%     
%     Inh_bin{sess}=inh_bin;
%     Inh_bin2{sess}=inh_bin2;
%     
%     FT{sess}=frametrigger;
%     Inh_onset{sess}=inh_onset;
%     FvOnTime{sess}=fvOnTime;
% end
% Sniff=cell2mat(sniff);
% Sniff_time=cell2mat(sniff_time);

for sess=sess_used
[sniff{sess,1},~,Data{sess},~]=read_sniff_frametrigger_trialinfo(H5{sess},path_h5,pre,post);
end
Sniff=cell2mat(sniff(cellfun(@(x) ~isempty(x),sniff)));
%%
odor_all=repmat(odor_num,1,numel(H5));%concatenate odor_num for all sessiosn

% if length(frametrigger)<size(SV{1},1)
%    SV{1}(1:size(SV{1},1)-length(frametrigger),:)=[];
% end
%%
%Check sniff alignment to fvOnTime or inh_onset
%{
for sess=sess_used
inh_bin=findfirst(sniff_time{sess}==double(inh_onset),2);
fv_bin=findfirst(sniff_time{sess}==double(fvOnTime),2);


offset=mean(sniff{sess}(:));
pn=sniff{sess}>offset;
filt=[ones(1,30),-1*ones(1,30)];
trig=zeros(size(sniff{sess}));
for t=1:size(sniff{sess},2)-length(filt)
    for tr=1:size(sniff{sess},1)
        trig(tr,t+length(filt)/2)=pn(tr,t:t+length(filt)-1)*filt';
        trig(tr,1:fv_bin(tr)+10)=0;
    end
end
inh_bin2=findfirst(trig==30,2);

for tr=1:size(sniff{sess},1)
   sn_inh(:,tr)=sniff{sess}(tr,inh_bin(tr):inh_bin(tr)+300);
   sn_inh2(:,tr)=sniff{sess}(tr,inh_bin2(tr):inh_bin2(tr)+300);
   sn_fv(:,tr)=sniff{sess}(tr,fv_bin(tr):fv_bin(tr)+300);
end

figure2;
srange=1:45;
subplot(221);plot(sn_inh(:,srange));
subplot(222);plot(sn_fv(:,srange));
subplot(223);plot(sn_inh2(:,srange));
Inh_bin{sess}=inh_bin;
Inh_bin2{sess}=inh_bin2;
end
%}
%%
clear Dff
for sess=sess_used
data=Data{sess};
    data.fps=30;
    [Dff{sess},TI{sess},trials_read]=ReconstructionSVD_2p(U,SV{sess},cellMask_vec,data,pre,post);%Gcamp
    
end
odor_all=repmat(odor_num,1,length(sess_used));%concatenate odor_num for all sessiosn
dff_all=cat(3,Dff{:});

%PID signal
%Calculate df/f for the first 1sec, 2sec, ...
[data,path_used]=read_h5(h5_name_PID,path_h5);
[PID,~,~,PID_time]=read_sniff_frametrigger_trialinfo(h5_name_PID,path_used,pre,post,false);
odor_num2=odor_num(2:end-1);%In PID session, I forgot to add empty trials at the beginning and end of session

if field<=2
%plot mean df/f for each reference odor
figure2('All cells');
ylm=[-0.6,0.5];
setym = @(x) set(gca,'YTick',-0.2:0.2:0.4,'YTickLabel',-0.2:0.2:0.4,'YLim',[-0.5,0.4]);
setylim=setym;
setym2 = @(x) set(gca,'YTick',-0.2:0.2:0.4,'YTickLabel',-0.2:0.2:0.4,'YLim',[-0.2,0.4]);
setylim2=setym2;

r=[100,90,71,50,29,10,0];
for i=1:8
    subplot(2,8,i)
    plot(mean(dff_all(:,:,odor_all==i),3)');hold on
    box off
    set(gca,'XTick',500:500:1500,'XTickLabel',0:0.5:1)
    setylim();
    line([pre,pre],ylim,'LineStyle','--','Color','k');
    line([pre+1000,pre+1000],ylim,'LineStyle','--','Color','k');
    yy=ylim;
%     plot(mean(Sniff(odor_all==i,:))/1500-abs(yy(1))*0.5);hold on
%     plot(mean(PID(odor_num==i,:))/300-abs(yy(1))*0.7);hold on
    plot(0.4*mean(Sniff(odor_all==i,:))/max(Sniff(:))-abs(yy(1))*0.6,'k');hold on
    plot(0.4*mean(PID(odor_num2==i,:))/max(PID(:))-abs(yy(1))*0.9,'c');hold on
    if i==1
        xlabel('sec')
        ylabel('df/f')
        title('Empty vial','FontSize',8)
    else
    title(sprintf('%dEB/%dPPA',100-r(i-1),r(i-1)),'FontSize',8)
    end
    
    subplot(2,8,i+8)
    plot(mean(dff_all(:,:,odor_all==i),3)');hold on
    box off
    line([pre,pre],ylim);
    xlim([451,900]);
    set(gca,'XTick',500:200:900,'XTickLabel',0:200:400);
    setylim();

    line([pre,pre],ylim,'LineStyle','--','Color','k');
plot(0.4*mean(Sniff(odor_all==i,:))/max(Sniff(:))-abs(yy(1))*0.6,'k');hold on
    plot(0.4*mean(PID(odor_num2==i,:))/max(PID(:))-abs(yy(1))*0.9,'c');hold on
    if i==1
        xlabel('ms')
        ylabel('df/f')
    end
end
save_figs_2p(sprintf('%sAllSess_AllCells_EB_PPA',name_data),sn)

figure2;
morph2=[1,9:15];
for j=1:8
    subplot(2,8,j)
    i=morph2(j);
    plot(mean(dff_all(:,:,odor_all==i),3)');hold on
    box off
    set(gca,'XTick',500:500:1500,'XTickLabel',0:0.5:1)
    setylim();
    line([pre,pre],ylim,'LineStyle','--','Color','k');
    line([pre+1000,pre+1000],ylim,'LineStyle','--','Color','k');
    yy=ylim;
%     plot(mean(Sniff(odor_all==i,:))/1500-abs(yy(1))*0.5);hold on
%     plot(mean(PID(odor_num==i,:))/300-abs(yy(1))*0.7);hold on
    plot(0.4*mean(Sniff(odor_all==i,:))/max(Sniff(:))-abs(yy(1))*0.6,'k');hold on
    plot(0.4*mean(PID(odor_num2==i,:))/max(PID(:))-abs(yy(1))*0.9,'c');hold on
    if j==1
        xlabel('sec')
        ylabel('df/f')
        title('Empty vial','FontSize',8)
    else
    title(sprintf('%dBut/%dBz',100-r(j-1),r(j-1)),'FontSize',8)
    end
    
    subplot(2,8,j+8)
    plot(mean(dff_all(:,:,odor_all==i),3)');hold on
    box off
    line([pre,pre],ylim);
    xlim([451,900]);
    set(gca,'XTick',500:200:900,'XTickLabel',0:200:400);
    setylim();

    line([pre,pre],ylim,'LineStyle','--','Color','k');
plot(0.4*mean(Sniff(odor_all==i,:))/max(Sniff(:))-abs(yy(1))*0.6,'k');hold on
    plot(0.4*mean(PID(odor_num2==i,:))/max(PID(:))-abs(yy(1))*0.9,'c');hold on
    if j==1
        xlabel('ms')
        ylabel('df/f')
    end
end
save_figs_2p(sprintf('%sAllSess_AllCells_ButyricAcid_Benzaldehyde',name_data),sn)


%Separate plot for all odors/individual cells
clear dff_mean
for i=1:15
    dff_mean(:,:,i)=mean(dff_all(:,:,odor_all==i),3);
end

col=[0,0,0;winter(7);cool(7)];
op={[1:8],[1,9:15]};
odset={'EB-PPA','BzAld-ButA'};
for p=1:2
figure2(odset{p});
od_plot=op{p};
for i=1:min(num_cell,42)
    subplot(6,7,i)
    set(gca, 'ColorOrder', col(od_plot,:), 'NextPlot', 'replacechildren');
    
    plot(squeeze(dff_mean(i,:,od_plot)));hold on
    box off
    set(gca,'XTick',500:500:1500,'XTickLabel',{'0','0.5','1'})
    
    setylim2();
    line([pre,pre],ylim,'LineStyle','--','Color','k');
    line([pre+1000,pre+1000],ylim,'LineStyle','--','Color','k');
%     if i==31
%         xlabel('(sec)');ylabel('df/f')
%     end
end
hL=legend(odor_set(od_plot),'FontSize',6);
newPosition = [0.94 0.6 0.03 0.2];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
save_figs_2p(sprintf('%sAllOdors_IndCell_%s',name_data,odset{p}),sn)
end

col=[0,0,0;winter(7);cool(7)];
op={[1:8],[1,9:15]};
for p=1:2
figure2(odset{p});
od_plot=op{p};
for i=1:min(num_cell,42)
    subplot(6,7,i)
    set(gca, 'ColorOrder', col(od_plot,:), 'NextPlot', 'replacechildren');
    
    plot(squeeze(dff_mean(i,:,od_plot)));hold on
    box off
    setylim2();
    line([pre,pre],ylim);
    xlim([451,900]);
    set(gca,'XTick',500:200:900,'XTickLabel',0:200:400);
%     if i==31
%         xlabel('(ms)');ylabel('df/f')
%     end
end
legout(odor_set(od_plot))
save_figs_2p(sprintf('%sAllOdors_IndCell_AroundInh_%s',name_data,odset{p}),sn)

end



%Single trial analysis with sniff
col=[0,0,0;winter(7);cool(7)];
subt=@(m,n,p) subtightplot(m,n,p,[0.05,0.05],[0.03,0.03],[0.03,0.03]);
%  choose subset of cells sensitive to EB
maxf=abs(mean(max(dff_all(:,500:1500,odor_all==2|odor_all==8),[],2),3));
ind=maxf>0.15;
for sess=sess_used
    figure2(sprintf('sess%d',sess));
    for tr=1:64
        subt(8,8,tr)
        plot(Dff{sess}(ind,:,tr+1)','Color',col(odor_num(tr+1),:));hold on
        plot((sniff{sess}(tr,:))/5000-0.35);hold on
        xlim([451,950]);set(gca,'XTick',500:200:900,'XTickLabel',0:200:400)
        line([pre,pre],ylim)
        setylim();
        box off
        title(sprintf('Trial%d',tr),'FontSize',8)
    end
end
% legout(odor_set)
% save_figs_2p(sprintf('%sSubset_SingleTrialDff_Sniff',name_data),sn)

elseif field==3
    odor_num=[1,1,2*ones(1,15),1,3*ones(1,15),1];

    odor_all=odor_num;%concatenate odor_num for all sessiosn

%plot mean df/f for each reference odor

ylm=[-0.6,0.6];
setym = @(x) set(gca,'YTick',-0.2:0.2:0.6,'YTickLabel',-0.2:0.2:0.6,'YLim',[-0.5,0.6]);
setylim=setym;
setym2 = @(x) set(gca,'YTick',-0.2:0.2:0.6,'YTickLabel',-0.2:0.2:0.6,'YLim',[-0.2,0.6]);
setylim2=setym2;

%Single trial analysis with sniff
col=[0,0,0;0,0,1;1,0,0];
subt=@(m,n,p) subtightplot(m,n,p,[0.05,0.05],[0.03,0.03],[0.03,0.03]);
%  choose subset of cells sensitive to EB
maxf=abs(mean(max(dff_all(:,500:1500,odor_all==2|odor_all==8),[],2),3));
ind=maxf>0.2;
for sess=sess_used
    figure2(sprintf('sess%d',sess));
    for tr=1:32
        subt(6,6,tr)
        plot(Dff{sess}(ind,:,tr+1)','Color',col(odor_num(tr+1),:));hold on
        plot((sniff{sess}(tr,:))/5000-0.35);hold on
        xlim([451,950]);set(gca,'XTick',500:200:900,'XTickLabel',0:200:400)
        line([pre,pre],ylim)
        setylim();
        box off
        title(sprintf('Trial%d',tr),'FontSize',8)
    end
end

%All trials/cell
figure2;
col=[0,0,0;winter(15)];
for i=1:num_cell
    subt(6,7,i)
    set(gca, 'ColorOrder', col, 'NextPlot', 'replacechildren');
    plot(squeeze(Dff{sess}(i,:,18:33)));hold on
        plot((sniff{sess}(tr,:))/5000-0.35);hold on
        xlim([451,950]);set(gca,'XTick',500:200:900,'XTickLabel',0:200:400)
        line([pre,pre],ylim)
        setylim();
        box off

end

end



%%
% pca

res_amp=max(max(abs(dff_mean(2:end,501:1000,:)),[],3),[],2);
cell_used = find(res_amp>prctile(res_amp,50));
dma = dff_mean(cell_used,:,:);
[coeff,score,latent]=pca(reshape(dma,size(dma,1),[])');

sall = reshape(score(:,1:3)',[],size(dma,2),size(dma,3));
col=[0,0,0;summer(7);cool(7)];
figure2;
trange=401:1000;
% ms = [2*ones(1,50),12*ones(1,1000),2*ones(1,500)];
ms = [2*ones(1,100),9*ones(1,500)];
for i=9:15
    
    scatter3(sall(1,trange,i),sall(2,trange,i),sall(3,trange,i),ms,col(i,:),'filled')
    hold on
end
legend(odor_set)


