function [frlips frvideowithlips frvideo dim ampix mov2] = getthelipstuff(vid, varargin)
addpath(genpath('F:\Other\Programs\Matlab\UGM'));

% first decide what the dimensions are
fr = 1;
if isempty(varargin)
    dat = rgb2gray(vid.frames(fr).cdata);
    [dattemp dim] = imcrop(dat);
    dim = [round(dim(2)) round(dim(2)+dim(4)); round(dim(1)) round(dim(1)+dim(3))];
    imagesc(dat(dim(1,1):dim(1,2), dim(2,1):dim(2,2)));
    pause
else
    dim = varargin{1};
end;

for fr = 1:vid.nrFramesTotal
    dat = rgb2gray(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),:));
    
    frvideo(fr,:,:) = dat;
    frvideowithlips(fr,:,:) = dat;
    
    H = 256.*(double(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),2))./double(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),1)));
    I = (double(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),1))+double(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),2))+double(vid.frames(fr).cdata(dim(1,1):dim(1,2), dim(2,1):dim(2,2),3)))./3;
    
    % Decoding with UGM_DECODE_ICMrestart
    UGM = gotoUGM(H);
    
    close all
    nRestarts = 500;
    ICMrestartDecoding = UGM_Decode_ICMrestart(UGM.nodePot,UGM.edgePot,UGM.edgeStruct,nRestarts);
    %ICMrestartDecoding = UGM_Decode_ICM(UGM.nodePot, UGM.edgePot, UGM.edgeStruct);
    
%     figure;
%     subplot(2,1,1)
%     imagesc(reshape(ICMrestartDecoding,UGM.nRows,UGM.nCols));
%     colormap gray
%     subplot(2,1,2)
%     imagesc(H);
%     title('ICM with Restarts Decoding of Noisy X');
    
    % % double check whether it found anything on the edges (this shouldn't
    % % happen) DOESNT WORK YET...
    % [xt r] = getcontours2(1:UGM.nCols, 1:UGM.nRows, double(reshape(ICMrestartDecoding-1, UGM.nRows, UGM.nCols)));
    % imagesc(xt)
    % colormap gray
    
    clear clus cntrs
    % only keep the biggest edge in there (the rest is noise)
    freqnm = [1:UGM.nRows];
    timenm = [1:UGM.nCols];
    h = figure;
    [xt] = contour(1:UGM.nCols, 1:UGM.nRows, squeeze(double(reshape(ICMrestartDecoding-1, UGM.nRows, UGM.nCols))),1, 'linecolor', [0 0 0]);
    close(h);
    locs = find(xt(1,:) == 0.5);
    for t = 1:length(locs)
        if t ~= length(locs)
            cntrs{t} = xt(:,locs(t)+1:locs(t+1)-1);
        else
            cntrs{t} = xt(:,locs(t)+1:end);
        end;
        clus.minx(t) = min(cntrs{t}(1,:));
        clus.maxx(t) = max(cntrs{t}(1,:));
        clus.miny(t) = min(cntrs{t}(2,:));
        clus.maxy(t) = max(cntrs{t}(2,:));
        clus.meanx(t) = mean(cntrs{t}(1,:));
        clus.meany(t) = mean(cntrs{t}(2,:));
        clus.difmean(t) = abs((clus.meanx(t)-UGM.nCols/2)) + abs((clus.meany(t)-UGM.nRows/2));
    end;
    %%
    newmask = squeeze(double(reshape(ICMrestartDecoding-1, UGM.nRows, UGM.nCols)));
    for t = 1:length(clus.minx)
        if clus.maxx(t)-clus.minx(t) < 20
            clus.difmean(t) = 100;            
        end;
        [X Y] = meshgrid(1:UGM.nCols,1:UGM.nRows);
        x = reshape(X,1,numel(X));
        y = reshape(Y,1,numel(Y));
        IN = inpolygon(x, y, cntrs{t}(1,:), cntrs{t}(2,:));
        clus.ampix(t) = length(find(IN == 1));
    end;
    [themax locmax] = max(clus.ampix);
    rmtrm = [];
    for t = 1:length(clus.minx)
        if (clus.ampix(t) > themax-40 & clus.ampix(t) < themax+40) & t ~= locmax 
            rmtrm = t;            
        end;
    end;
    ampix(fr) = themax;
    clus.ampix(rmtrm) = [];
    clus.difmean(rmtrm) =[];
    cntrs(rmtrm) = [];
        
    [temp ord] = sort(clus.ampix);
    ord = ord(end:-1:1);
    %[temp lipscont] = max(clus.ampix);
    
    for t = 1:length(cntrs)
        if t ~= 1;
            %newmask(floor(clus.minx(t)):ceil(clus.maxx(t)),floor(clus.miny(t)):ceil(clus.maxy(t))) = 0;
            %newmask(floor(clus.miny(t)):ceil(clus.maxy(t)),floor(clus.minx(t)):ceil(clus.maxx(t))) = 1;
            [X Y] = meshgrid(1:UGM.nCols,1:UGM.nRows);
            x = reshape(X,1,numel(X));
            y = reshape(Y,1,numel(Y));
            IN = inpolygon(x, y, cntrs{ord(t)}(1,:), cntrs{ord(t)}(2,:));
            subplot(2,1,1)
            hold on
            plot(x(IN),y(IN),'r+');
            tempINmatrix = squeeze(double(reshape(IN,  UGM.nRows, UGM.nCols)));            
            newmask(tempINmatrix == 1) = 1;
        else
            [X Y] = meshgrid(1:UGM.nCols,1:UGM.nRows);
            x = reshape(X,1,numel(X));
            y = reshape(Y,1,numel(Y));
            IN = inpolygon(x, y, cntrs{ord(t)}(1,:), cntrs{ord(t)}(2,:));
            subplot(2,1,1)
            hold on
            plot(x(IN),y(IN),'b+');
            tempINmatrix = squeeze(double(reshape(IN,  UGM.nRows, UGM.nCols))); 
            newmask( tempINmatrix == 1) = 0;
        end;       
    end;
    imagesc(newmask)
    colormap gray
    
    %%
    dattemp = squeeze(frvideowithlips(fr,:,:));
    dattemp(find(newmask == 0)) = 1;
    frvideowithlips(fr,:,:) = dattemp;
    frlips(fr,:,:) = ICMrestartDecoding;
    
    mov(fr).cdata(:,:,1) = dattemp;
    mov(fr).cdata(:,:,2) = dattemp;
    mov(fr).cdata(:,:,3) = dattemp;
    mov(fr).colormap = [];
    display(['finished with ' num2str(fr) ' of ' vid.nrFramesTotal]);
end;
mov2 = mov;
resizeval = 5;
newCdata = cellfun(@(x) x(...
    repmat(1:size(x,1),resizeval,1), ...         % like kron, but then
    repmat(1:size(x,2), resizeval,1), :), ...    % a bit faster, and suited
    {mov.cdata}, 'UniformOutput', false);  % for 3D arrays

% assign all new data back to the movie
[mov2.cdata] = newCdata{:};

close all
hf = figure;
set(hf, 'position', [400 200 600 400])
movie(hf, mov2, 1, vid.rate);




