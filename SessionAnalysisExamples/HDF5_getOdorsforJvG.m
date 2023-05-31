function [odorInfo] = HDF5_getOdorsforJvG(fpathH5,fnameH5,trials_read, latfactor, upper,sortflow)

if ~exist('latfactor','var')
    latfactor= 10;
end
if ~exist('sortflow','var')
    sortflow=false;
end
if ~exist('upper','var')
    upper=true;
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
    dilution_factor = 1-round(H5.dilutors0x3Adilutor_00x3Avac_flow  ./10)./100;
end
rackOdors0 = deblank(string((H5.olfas0x3Aolfa_00x3Aodor)'));
rackOdors1 = deblank(string((H5.olfas0x3Aolfa_10x3Aodor)'));
if sortflow
    olfa0_flow = H5.olfas0x3Aolfa_00x3Amfc_1_flow;
    olfa1_flow = H5.olfas0x3Aolfa_10x3Amfc_1_flow;
else
    olfa0_flow = H5.olfas0x3Aolfa_00x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_00x3Avialconc;
    olfa1_flow = H5.olfas0x3Aolfa_10x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_10x3Avialconc;
end

if upper
    olfa0_flow(rackOdors0 =='None') = [];
    olfa1_flow(rackOdors1 =='None') = [];
    rackOdors0(rackOdors0 =='None') = '';
    rackOdors1(rackOdors1 =='None') = '';
else
    olfa1_flow(rackOdors0 == "empty" & rackOdors1 == "None") = olfa0_flow(rackOdors0 == "empty");
    rackOdors1(rackOdors0 == "empty" & rackOdors1 == "None") = rackOdors0(rackOdors0 == "empty");
    rackOdors0(rackOdors0 == "empty") =repmat('None', sum(rackOdors0 == "empty"),1);
    olfa0_flow(rackOdors0 =="empty") = [];
    olfa1_flow(rackOdors1 =='None') = [];
    olfa0_flow(rackOdors0 =="None") = [];
    rackOdors0(rackOdors0 =='None') = '';
    rackOdors1(rackOdors1 =='None') = '';
    rackOdors0(rackOdors0 == "empty") = '';
end

if sortflow
    olfa_odor_ = strcat(rackOdors0,rackOdors1, ':', num2str(olfa0_flow,'%.3f'),num2str(olfa1_flow,'%.3f'),' flow', StimID) ;

else
    olfa_odor_ = strcat(rackOdors0,rackOdors1, ':', num2str(olfa0_flow,'%.3f'),num2str(olfa1_flow,'%.3f'),' nM', StimID) ;
end

olfa_odors = olfa_odor_(trials_read);
odors = unique(olfa_odors);
for idx = 1:length(odors)
    odorTrials{idx} = find(olfa_odors==odors(idx));
end

odorInfo.odors = odors;
odorInfo.odorTrials = odorTrials;
end