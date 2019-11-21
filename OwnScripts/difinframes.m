function x = difinframes(stimvid)

[vid aud] = mmread(stimvid);
for i = 1:length(vid.frames)
    vidFrames(:,:,:,i) = vid.frames(i).cdata;
end;
numFrames = vid.nrFramesTotal;
FrameRate = vid.rate;
vidtimes = vid.times;

dims = [1 288;1 360]; 

% for every frame calculate the absolute difference between
% adjacent frames
for frms = 2:numFrames
    abdif(:,:,:,frms) = imabsdiff(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms-1), vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms));
    t = squeeze(abdif(:,:,:,frms));
    t = double(mean(t,3));
    t = t(:);
    %t = double(t(:,:,2));
    vardif(frms) = (var(t));
    meandif(frms) = mean(mean(t));
    mov(frms-1).cdata = vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms);
    mov(frms-1).colormap = [];
end;
%      hf = figure;
%      set(hf, 'position', [600 600 300 300])
%      % Playback movie once at the video's frame rate
%      movie(hf, mov, 1, FrameRate);
%       peaks.aud(am) = findpeaks(aud(am), 'Minpeakheight', 0.1, 'MINPEAKDISTANCE', 100000);
%      


[t loc] = max(aud(am,:));
peaks.aud{i}(am) = audtime(loc);
[t loc] = max(vardif{i}(am,:));
peaks.vis{i}(am) = vidtime(loc);
difs(end+1) = peaks.aud{i}(am)-peaks.vis{i}(am);

cntplot = cntplot + 1;
subplot(6,3,cntplot)
plot(vardif{i}(am,:));
title([nm{i} ' ' num2str(am)]);


vardifavg{i} = mean(vardif{i},1);
meandifavg{i} = mean(meandif{i},1);
audavg{i} = mean(aud{i},1);
