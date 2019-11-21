cd C:\Users\Sanne.tenOever\Documents\Experiments\AVspeechthings\
vidloc = 'C:\Users\Sanne.tenOever\Documents\Experiments\AVspeechthings\GRID\s2\video\mpg_6000\';
addpath('C:\Users\Sanne.tenOever\Documents\Experiments\AVspeechthings\mmread');
addpath('C:\Users\Sanne.tenOever\Documents\Experiments\AVspeechthings\detectfaceparts');

%%
vidname = 'bbaf1n.mpg';
vidstim = [vidloc vidname];

[vid aud] = mmread(vidstim);

[frlips frvideowithlips frvideo dim ampix mov] = getthelipstuff(vid);
close all
hf = figure;
set(hf, 'position', [400 100 800 800])
movie(hf, mov, 1, vid.rate);

subplot(2,1,1)
plot(ampix)
title('area of the mouth opening');
subplot(2,1,2)
plot(smooth(diff(ampix),5), 'r')
title('velocity');

%% check the frequency of visual signal
addpath('C:\Users\Sanne.tenOever\Documents\Programs\Matlab\fieldtrip-20120904\fieldtrip-20131208');
clear dat
dat.dimord = 'chan_time';
dat.trial{1}(1,:) = diff(ampix);
dat.label{1} = '1';
dat.time{1} = [0:vid.rate:vid.rate*73]./1000;

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
fftvid = ft_freqanalysis(cfg, dat);

figure
ft_singleplotER([], fftvid);

%% fft of aud
clear dat
dat.dimord = 'chan_time';
dat.trial{1}(1,:) = aud.data(:,1);
dat.label{1} = '1';
dat.time{1} = [0:1/44100:(length(aud.data)-1)/44100];

cfg =[];
cfg.detrend = 'no';
cfg.resamplefs = vid.rate;
dat = ft_resampledata(cfg, dat); 


