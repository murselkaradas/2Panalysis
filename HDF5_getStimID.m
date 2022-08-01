function [odorInfo] = HDF5_getStimID(fpathH5,fnameH5,trials_read, latfactor, sortflow)

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

PatternID = deblank(string((H5.patternid)'));
OffDur = H5.pulseOffDur_1;

PatternID(PatternID == "") = 'Blank';
%OffDur(strcmp(PatternID, "NoOverlap.tif")) = [];
%OffDur(strcmp(PatternID, "NoOverlapReversed.tif")) = [];
PatternID = PatternID +"--"+ num2str(OffDur)
PatternIDs = PatternID(trials_read);
IDs = unique(PatternIDs);
for idx = 1:length(IDs)
    IDsTrials{idx} = find(PatternIDs==IDs(idx));
end

odorInfo.odors = IDs;
odorInfo.odorTrials = IDsTrials;
end