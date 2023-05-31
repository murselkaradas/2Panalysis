function [odorInfo] = HDF5_getOdors(fpathH5,fnameH5,trials_read, latfactor, sortflow)

if ~exist('latfactor','var')
    latfactor= 10;
end
if ~exist('sortflow','var')
    sortflow=false;
end
H5=h5read(fullfile(fpathH5,fnameH5),'/Trials');

Ntrial = length(H5.trialNumber);
%Check trials read exist or not
if ~exist('trials_read','var')
    trials_read = logical(ones(1,Ntrial));
end

if any(ismember(fields(H5), 'amplitude_1'))
    StimID = strings([Ntrial,1]);
    stimon = find(H5.amplitude_1 >0);
%     lat = round((H5.laserontime(stimon) - H5.inh_onset(stimon))/latfactor)*latfactor;
    lat = round(H5.pulseOnsetDelay_1(stimon)/latfactor)*latfactor;
    StimID(stimon) = strcat('-SL', num2str(lat, '%d'));
    
end    
    
if ~any(ismember(fields(H5),'dilutors0x3Adilutor_00x3Aair_flow'))
    dilution_factor = ones(Ntrial,1);
else 
    dilution_factor = 1-round(H5.dilutors0x3Adilutor_00x3Aair_flow./10)./100;
end
rackOdors0 = deblank(string((H5.olfas0x3Aolfa_00x3Aodor)'));
rackOdors1 = deblank(string((H5.olfas0x3Aolfa_10x3Aodor)'));
if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
    rackOdors2 = deblank(string((H5.olfas0x3Aolfa_20x3Aodor)'));
end
if sortflow
    olfa0_flow = H5.olfas0x3Aolfa_00x3Amfc_1_flow;
    olfa1_flow = H5.olfas0x3Aolfa_10x3Amfc_1_flow;
    if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
        olfa2_flow = H5.olfas0x3Aolfa_20x3Amfc_0_flow;
    end
else
    olfa0_flow = H5.olfas0x3Aolfa_00x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_00x3Avialconc;
    olfa1_flow = H5.olfas0x3Aolfa_10x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_10x3Avialconc;
    if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
        olfa2_flow = H5.olfas0x3Aolfa_20x3Amfc_0_flow.*dilution_factor.*H5.olfas0x3Aolfa_20x3Avialconc;
    end
end

olfa0_flow(rackOdors0 =='None') = 0;
olfa1_flow(rackOdors1 =='None') = 0;
rackOdors0(rackOdors0 =='None') = '--';
rackOdors1(rackOdors1 =='None') = 'xx';
if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
    olfa2_flow(rackOdors2 =='None') = 0;
    rackOdors2(rackOdors2 =='None') = 'zz';
end
if sortflow
    if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
        olfa_odor_ = strcat(rackOdors0,'/',rackOdors1, '/',rackOdors2, ':', num2str(olfa0_flow,'%.3f'),'/',num2str(olfa1_flow,'%.3f'),'/',num2str(olfa2_flow,'%.3f'),'/flow','/', StimID);
    else
        olfa_odor_ = strcat(rackOdors0,'/',rackOdors1, ':', num2str(olfa0_flow,'%.3f'),'/',num2str(olfa1_flow,'%.3f'),'/flow','/', StimID) ;
    end
else
    if any(ismember(fields(H5),'olfas0x3Aolfa_20x3Aodor'))
        olfa_odor_ = strcat(rackOdors0,'/',rackOdors1,'/',rackOdors2, ':', num2str(olfa0_flow,'%.3f'),'/',num2str(olfa1_flow,'%.3f'),'/', num2str(olfa2_flow,'%.3f'),'/nM', '/',StimID) ;

    else
        olfa_odor_ = strcat(rackOdors0,'/',rackOdors1, ':', num2str(olfa0_flow,'%.3f'),'/',num2str(olfa1_flow,'%.3f'),'/nM', '/',StimID) ;

    end
end

olfa_odors = olfa_odor_(trials_read);
odors = unique(olfa_odors);
for idx = 1:length(odors)
    odorTrials{idx} = find(olfa_odors==odors(idx));
end

odorInfo.odors = odors;
odorInfo.odorTrials = odorTrials;
end
