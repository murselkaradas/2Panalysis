function stimframe_ =    findmostrepeatedstims(stimframe)
    stimframe_  = stimframe;
    for i = 1: size(stimframe,2)
        x =stimframe(:,i);
        [n,bin] = hist(x,unique(x));
        [~, idx] = sort(-n);
        stimtimes= bin(idx(1:2));
        if abs(stimtimes(1) - stimtimes(2)) ==1
            stimframe_((stimframe(:,i) ~= stimtimes(1) & stimframe(:,i) ~= stimtimes(2)), i) = max(stimtimes);
        
        else
           stimframe_(stimframe(:,i) ~= stimtimes(1),i) = stimtimes(1);
        end

    end

end