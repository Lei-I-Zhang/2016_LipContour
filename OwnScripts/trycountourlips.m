function x = trycountourlips(stimvid)
addpath(genpath('F:\Experiments\AVspeechthings'));
stimvid = 'F:\Experiments\AVspeechthings\s1_mediumquality\Video\bbaf3s.mpg';
[vid aud] = mmread(stimvid);
for i = 1:length(vid.frames)
    vidFrames(:,:,:,i) = vid.frames(i).cdata;
end;
numFrames = vid.nrFramesTotal;
FrameRate = vid.rate;
vidtimes = vid.times;

dims = [1 340;1 360];
dims = [200 250; 100 210];

image(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),:,1));

% (estimated) logarithmic color transform:
% first only for the first frame to see the hue of the lips.
H = 256.*(double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,1))./double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,1)));
I = (vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,1)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,1)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),3,1))./3;

subplot(3,1,1)
imagesc(rgb2gray(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),:,1)));
subplot(3,1,2);
imagesc(H)
set(gca, 'clim', [0 300]);
colormap('gray');
subplot(3,1,3);
imagesc(I)
set(gca, 'clim', [0 200]);
colormap('gray');

% extract a lip point:
Xvalue = 62;
Yvalue = 26;
lipH = H(Yvalue, Xvalue);

%% estimate Hlip deltaH and thetaH
hist(H(:), 200);
[modeH] = mode(H(:));
T = H >=modeH-3.5 & H <=modeH+3.5;
for yax = 1:size(H,1) % the different y pixels
    for xax = 1:size(H,2) % the different x pixels
        if abs(H(yax, xax)-modeH)./std(H(T)) <= 16 % !!!not sure whether the variance here is based on the single pixel or whole data
            Hnew(yax, xax) = (256-(H(yax, xax)-modeH)./std(H(T)))^2; % filter to value of lip
        else
            Hnew(yax, xax) = 0;
        end;
    end;
end;

subplot(2,1,1)
imagesc(rgb2gray(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),:,1)));
subplot(2,1,2)
imagesc(Hnew)
colormap('gray');
colorbar

%% same for the lips now
T = Hnew ~= 0;
Hnoskin = H;
Hnoskin(T) = NaN;
imagesc(Hnoskin);
colormap('gray');
colorbar
hist(Hnoskin(:), 50);
T = Hnoskin > 300;
Hnoskintemp = Hnoskin;
Hnoskintemp(T) = NaN;

imagesc(Hnoskintemp)
colormap('gray');
colorbar

hist(Hnoskintemp(:), 50);
[modeHlip] = mode(Hnoskintemp(:));
deltaHlip = 8;

%% define the threshold to suppress camera noise
for frms = 2:numFrames
    I = double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,frms)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,frms)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),3,frms))./3;
    Imin1 = double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,frms-1)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,frms-1)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),3,frms-1))./3;
    for yax = 1:size(H,1)
        for xax = 1:size(H,2)
            fd(frms-1,yax,xax) = abs(I(yax,xax)-Imin1(yax,xax));
        end;
    end;
    theta(frms-1) = 2^(entropy(fd(frms-1,:,:)));
end;

%% try to do the thingy....
% label motion and red hueness...
for frms = 2:numFrames
    H = 256.*(double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,frms))./double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,frms)));
    for yax = 1:size(H,1)
        for xax =1:size(H,2)
            if fd(frms-1,yax,xax) > theta(frms-1)
                motion(frms-1,yax,xax) = 1;
            else
                motion(frms-1,yax,xax) = 0;
            end;
            if H(yax,xax) > modeHlip - deltaHlip && H(yax,xax) < modeHlip + deltaHlip
                hue(frms-1,yax,xax) = 1;
            else
                hue(frms-1,yax,xax) = 0;
            end;
            if motion(frms-1,yax,xax)== 0 && hue(frms-1,yax,xax) == 0
                label(frms-1, yax, xax) = 1;
            elseif motion(frms-1,yax,xax)== 1 && hue(frms-1,yax,xax) == 0
                label(frms-1, yax, xax) = 2;
            elseif motion(frms-1,yax,xax)== 0 && hue(frms-1,yax,xax) == 1
                label(frms-1, yax, xax) = 3;
            elseif motion(frms-1,yax,xax)== 1 && hue(frms-1,yax,xax) == 1
                label(frms-1, yax, xax) = 4;
            end;
        end;
    end;
end;

%% maximizing a posteriori probablity of the label field (minimizing a global energy function:
% so itterative process minimizing the following:
% neighbours:
values = [0 0 0 1 1 1 -1 -1 -1];
neighbours = [];
while size(neighbours,1) ~= 26
    i = size(neighbours,1)+1;
    neighbours(i,:) = values(random('unid', [9 9 9]));
    c = 0;
    for ij = 1:size(neighbours,1)-1
        if ((neighbours(ij,1) == neighbours(end,1) && neighbours(ij,2) == neighbours(end,2)) && neighbours(ij,3) == neighbours(end,3)) && c == 0
            neighbours(end,:) = [];
            c = 1;
        end;
    end;
    %display(size(t,1));
end;

lamba = 1;
alpha = 20;

for it = 1:20
    for frms = 2:numFrames-2
        for yax = 1:size(H,1)
            for xax =1:size(H,2)
                for n = 1:size(neighbours,1)  
                    Bx = (1+abs(hue(frms,xax,yax)-hue(frms+neighbours(n,1),xax,yax)) + abs(motion(frms,xax,yax)-hue(frms+neighbours(n,1),xax,yax)))^2;
                    By = ((1+abs(hue(frms,xax,yax)-hue(frms,xax+neighbours(n,2),yax)) + abs(motion(frms,xax,yax)-hue(frms,xax+neighbours(n,2),yax)))*2)^2;
                    Btcal = sum(hue(frms,xax,yax), hue(frms+neighbours(n,3), xax, yax), motion(frms,xax,yax));
                    if Btcal == 0 || Btcal == 2
                        Bt = 1;
                    elseif Btcal == 1 || Btcal == 3
                        Bt = 2;
                    end;                    
                    V(yax,xax,n) = (Bx*Bt)./(sqrt(Bt^2*(neighbours(2)^2+4*neighbours(3)^2)+Bs^2*neighbours(1)^2));
                end;
                V2(yax,xax) = sum(V,3);
            end;
        end;
        Um(frms-1) = sum(V2(:));
    end;
    
    
    
    Vst = 1./sqrt();
    W(S) = (lamba*Uh(S)+Ufd(S))+ alpha*Um(S);
end;

W(s) = sum(
U0(s) = sum(abs(abs(o


%%


T = H >=modeH-3 & H <=modeH+3;
for yax = 1:size(H,1) % the different y pixels
    for xax = 1:size(H,2) % the different x pixels
        if abs(Hnoskin(yax, xax)-modeHnoskin)./std(Hnoskin(T)) <= 16 % !!!not sure whether the variance here is based on the single pixel or whole data
            Hnew(yax, xax) = (256-(Hnoskin(yax, xax)-modeHnoskin)./std(Hnoskin(T)))^2; % filter to value of lip
        else
            Hnew(yax, xax) = 0;
        end;
    end;
end;

subplot(2,1,1)
imagesc(rgb2gray(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),:,1)));
subplot(2,1,2)
imagesc(Hnew)
colormap('gray');
colorbar
%%


H = squeeze(256.*(double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,:))./double(vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,:))));
I = double(squeeze((vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),1,:)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),2,:)+vidFrames(dims(1,1):dims(1,2), dims(2,1):dims(2,2),3,:))./3));

for neighbours = 1:size(H,3) % the different time points
    for yax = 1:size(H,1) % the different y pixels
        for xax = 1:size(H,2) % the different x pixels
            if abs(H(yax, xax,neighbours)-lipH)./std(H(:,:,neighbours)) <= 16 % !!!not sure whether the variance here is based on the single pixel or whole data
                Hnew(yax, xax,neighbours) = 0;
            else
                Hnew(yax, xax,neighbours) = (256-(H(yax, xax,neighbours)-lipH)./std(H(:,:,neighbours)))^2; % filter to value of lip
            end;
            if neighbours > 1
                fd(yax,xax,neighbours-1) = abs(I(yax,xax,neighbours)-I(yax,xax,neighbours-1)); % difference in intensity over time
            end
        end;
    end;
end;


