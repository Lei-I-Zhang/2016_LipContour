nm = {'Pa'; 'Ba'; 'Ta';'Da'; 'Ka';'Ga'};
dims = [190,230;100,180];
difs = [];
cntplot = 0;

cam = 2; %which it is for cut
for i = 1:6 % the stimuli
    for am = 1:3
        locstimuli = 'E:\Experiments\VowelExp\Own-Made-vowel-stimuli\endstim\Top_6_perstimulus\3chosen\';
        readerobj = VideoReader([locstimuli 'KH' nm{i} '_' num2str(am) '_fade.avi']);
        vidFrames = read(readerobj);
        numFrames = get(readerobj, 'NumberOfFrames');
        FrameRate = get(readerobj, 'FrameRate');
        
        [aud{i}(am,:),fs, nbits] = wavread([locstimuli 'KH' nm{i} '_' num2str(am) '_600ms.wav']);
        
        vidtime = [0:1/FrameRate:readerobj.Duration];
        vidtime(end) = [];
        
        audtime = [0:1/fs:.6];
        audtime = audtime + (1/FrameRate)*44-.4; % at 43th frame is the max aud in visual and .4 in the aud (however because other is difference, one frame later used
        
        % for every frame calculate the absolute difference between
        % adjacent frames
        for frms = 2:numFrames
            abdif{i}(am,:,:,:,frms) = imabsdiff(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms-1), vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms));
            t = squeeze(abdif{i}(am,:,:,:,frms));
            t = double(mean(t,3));
            t = t(:);
            %t = double(t(:,:,2));
            vardif{i}(am,frms) = (var(t));
            meandif{i}(am,frms) = mean(mean(t));
            mov(frms-1).cdata = vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2), :, frms);
            mov(frms-1).colormap = [];
        end;
        %      hf = figure;
        %      set(hf, 'position', [600 600 300 300])
        %      % Playback movie once at the video's frame rate
        %      movie(hf, mov, 1, readerobj.FrameRate);     
        %       peaks.aud{i}(am) = findpeaks(aud{i}(am), 'Minpeakheight', 0.1, 'MINPEAKDISTANCE', 100000); 
        [t loc] = max(aud{i}(am,:)); 
        peaks.aud{i}(am) = audtime(loc);
        [t loc] = max(vardif{i}(am,:));
        peaks.vis{i}(am) = vidtime(loc);
        difs(end+1) = peaks.aud{i}(am)-peaks.vis{i}(am);
        
        cntplot = cntplot + 1;
        subplot(6,3,cntplot)
        plot(vardif{i}(am,:));
        title([nm{i} ' ' num2str(am)]);
    end;
    vardifavg{i} = mean(vardif{i},1);
    meandifavg{i} = mean(meandif{i},1);
    audavg{i} = mean(aud{i},1);     
    
    %figure
    %plot(vardifavg{i});
end;

%% plot
figure
for i = 1:6
    subplot(3,2,i)
    plot(audtime, audavg{i}*50, 'r')
    hold on
    plot(vidtime(16:61), vardifavg{i}(16:61), 'Linewidth', 2)
    title(nm{i}, 'Fontsize', 15)
    set(gca, 'Fontsize', 15);
end;

%% average of the location of articulation
aloc = {'bilabial', 'coronal', 'dorsal'};
for i = 1:3
    vardifavg2{i} = mean([vardif{i*2-1}; vardif{i*2}]);
    meandifavg2{i} = mean([meandif{i*2-1}; meandif{i*2}]);
    audavg2{i} = mean([aud{i*2-1}; aud{i*2}]);
    subplot(3,1,i)
    plot(audtime, audavg2{i}*50, 'r')
    hold on
    plot(vidtime(16:61), vardifavg2{i}(16:61), 'Linewidth', 2)
    title(aloc{i}, 'Fontsize', 15)
    set(gca, 'Fontsize', 15);
end;


%%
addpath(genpath('C:\Users\Sanne.tenOever\Documents\Programs\Matlab\fieldtrip-20120904\fieldtrip-20120904'));
stim = 6;
rep = 1;
audtimeN = [0:1/fs:(audtime(1)-1/fs) audtime];
audsign = [zeros(1,length(0:1/fs:(audtime(1)-1/fs))) aud{stim}(rep,:)];
vissign = vardif{stim}(rep,:);

locstimuli = 'E:\Experiments\VowelExp\Own-Made-vowel-stimuli\endstim\Top_6_perstimulus\3chosen\';
readerobj = VideoReader([locstimuli 'KH' nm{stim} '_' num2str(rep) '_fade.avi']);
vidFrames = read(readerobj);
numFrames = get(readerobj, 'NumberOfFrames');
FrameRate = get(readerobj, 'FrameRate');

%%
close all
figure
set(gcf, 'position', [100 100 1000 200]);
d = 20;
for r = 1:5
    subplot(6,5,r)    
    imagesc(vidFrames(190:240,110:170,:,r+d));
    set(gca, 'xtick', 1000, 'ytick', 1000);
end;
subplot(10,1,3)
plot(vidtime(r+d:r+d+5),vissign(r+d:r+d+5));
title('velocity mouth movements');
set(gca,'xlim',[vidtime(r+d) vidtime(r+d+5)], 'xticklabel', '');
subplot(10,1,5)
t1 = nearest(audtimeN,vidtime(r+d));
t2 = nearest(audtimeN, vidtime(r+d+5));
plot(audtimeN(t1:t2), audsign(t1:t2));
title('auditory signal');
h = xlabel('time (sec)', 'Position', [vidtime(r+d+5) -3.2381 1.0001], 'HorizontalAlignment', 'right');
set(gca,'xlim',[vidtime(r+d) vidtime(r+d+5)]);
%exportfig(gcf, 'test2', 'color', 'rgb', 'preview', 'tiff');
loc
saveas(gcf, 'test.tif');



