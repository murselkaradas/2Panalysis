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
PulseDur = H5.pulseOnDur_1/1e3;
AOMpower = H5.AOMpower/1000;
PulseLat = H5.pulseOnsetDelay_1;
PatternID = deblank(string((H5.patternid)'));
OffDur = H5.pulseOffDur_1;

PulseDur(PatternID == "") = 0;
AOMpower(PatternID == "") = 0;
PulseLat(PatternID == "") = 0;
OffDur(PatternID == "") = 0;
PatternID(PatternID == "") = 'Blank';
cond1 = size(unique(PulseDur(PatternID ~= 'Blank')),1)>1;
cond2 = size(unique(OffDur(PatternID ~= 'Blank')),1)>1;
cond3 = size(unique(PulseLat(PatternID ~= 'Blank')),1)>1;
cond4 = size(unique(AOMpower(PatternID ~= 'Blank')),1)>1;
ID  = strings(Ntrial,1);
for i = 1:Ntrial
    if PatternID(i) =='Blank'
        ID(i) = strcat(ID(i) ,'id:',PatternID(i));
    else
        ID(i) = strcat(ID(i),'id:',PatternID(i));
        if cond1
            ID(i) = strcat(ID(i),'/PD:',num2str(PulseDur(i)));
        end
        if cond2
            ID(i) = strcat(ID(i),'/Off:',num2str(OffDur(i)));
        end
        if cond3
            ID(i) = strcat(ID(i),'/L:',num2str(PulseLat(i)));
        end
        if cond4
            ID(i) = strcat(ID(i),'/P:',num2str(AOMpower(i),'%.1f'));
        end
    end
end
%ID = strcat('ID:',PatternID,'/ONDur:',num2str(PulseDur),'/OFFDur:',num2str(OffDur),'/Lat:',num2str(PulseLat),'/Power:',num2str(AOMpower,'%.1f'));
IDs_ = ID(trials_read);
ID_unique = unique(IDs_);
for idx = 1:length(ID_unique)
    IDsTrials{idx} = find(IDs_==ID_unique(idx));
end

odorInfo.odors = ID_unique;
odorInfo.odorTrials = IDsTrials;
end