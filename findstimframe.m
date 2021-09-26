function stim_frames_updated =     findstimframe(data,stim_frames,halfwindow)
    % Check framebound is exist or not
    if ~exist('halfwindow','var')
         halfwindow = 1;
    end
   %stim_frames = Data.stim_frame(find(Data.stim_frame<= Framebound(2) & Data.stim_frame> Framebound(1)))-Framebound(1);
   
   if ~isempty(stim_frames)
       stim_frames_updated = stim_frames;
       for i = 1:length(stim_frames)
           frame_lbound = max(1, stim_frames(i)-halfwindow);
           frame_hbound = min(length(data), stim_frames(i)+halfwindow);
           [~,ind] = max(abs(data(frame_lbound:frame_hbound) - mean(data(frame_lbound-3:frame_hbound+3))));
           stim_frames_updated(i) = stim_frames(i) - halfwindow -1 + ind;

       end
   else
       stim_frames_updated = [];
   end


end

