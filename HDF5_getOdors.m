function [odorInfo] = HDF5_getOdors(fpathH5,fnameH5,trials_read)

H5=h5read(fullfile(fpathH5,fnameH5),'/Trials');

Ntrial = length(H5.trialNumber);
%Check trials read exist or not
if ~exist('trials_read','var')
    trials_read = logical(ones(1,Ntrial));
end

if any(ismember(fields(H5), 'amplitude_1'))
    StimID = strings([Ntrial,1]);
    stimon = find(H5.amplitude_1 >0);
    StimID(stimon) = '-S';
    
end    
    
if ~any(ismember(fields(H5),'dilutors0x3Adilutor_00x3Aair_flow'))
    dilution_factor = ones(Ntrial,1);
else 
    dilution_factor = (1-H5.dilutors0x3Adilutor_00x3Aair_flow./10000);
end
rackOdors0 = deblank(string((H5.olfas0x3Aolfa_00x3Aodor)'));
olfa0_flow = H5.olfas0x3Aolfa_00x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_00x3Avialconc;
rackOdors1 = deblank(string((H5.olfas0x3Aolfa_10x3Aodor)'));
olfa1_flow = H5.olfas0x3Aolfa_10x3Amfc_1_flow.*dilution_factor.*H5.olfas0x3Aolfa_10x3Avialconc;


olfa0_flow(rackOdors0 =='None') = 0;
olfa1_flow(rackOdors1 =='None') = 0;
rackOdors0(rackOdors0 =='None') = '---';
rackOdors1(rackOdors1 =='None') = '---';
olfa_odor_ = strcat(rackOdors0,'/',rackOdors1, ':', num2str(olfa0_flow,'%.2f'),'/',num2str(olfa1_flow,'%.2f'),' nM', StimID) ;

olfa_odors = olfa_odor_(trials_read);
odors = unique(olfa_odors);
for idx = 1:length(odors)
    odorTrials{idx} = find(olfa_odors==odors(idx));
end

odorInfo.odors = odors;
odorInfo.odorTrials = odorTrials;
end