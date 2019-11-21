function [tempmask clusideas] = getcontours2(x, y, mask, varargin);
freqnm = [1:size(y)];
timenm = [1:size(x)];
h = figure;
[xt] = contour(x, y, double(squeeze(mask)),1, 'linecolor', [1 1 1]);
close(h);
for t = 1:length(xt)-1 % check whether there are constant steps between the matrix dimensions
    tp(1,t)= abs(xt(1,t+1)-xt(1,t))<1;
    tp(2,t) = abs(xt(2,t+1)-xt(2,t))<1;
end;
allp = find(tp(1,:) == 0 & tp(2,:) == 0); %new edge as well as time or freq are not constant step
sizes = allp(2:end)- allp(1:end-1); % length of edge
sizes = [allp(1) sizes]; % first length is the number itself
sizes(2,:) = allp;
%abthres = find(sizes(1,:) > thres); % point that overgo the threshold.

tempmask = mask;
for z = 1:length(sizes(1,:))
    if z == 1
        examp = xt(:,1:sizes(2,z));
        clusideas.cenfreq(z) = mean(examp(1,:));
        clusideas.centime(z) = mean(examp(2,:));
    else
        examp = xt(:,sizes(2,z-1)+1:sizes(2,z));
        clusideas.cenfreq(z) = mean(examp(1,:));
        clusideas.centime(z) = mean(examp(2,:));
    end;
    mint = find(min(examp(1,:)) < x, 1, 'first');
    maxt = find(max(examp(1,:)) > x,1,'last');
    minf = find(min(examp(2,:)) < y, 1, 'first');
    maxf = find(max(examp(2,:)) > y, 1, 'last');
    if maxf > length(y)
        maxf = 79;
    end;
    if (x(maxt)-x(mint)) < 20
        tempmask(minf:maxf, mint:maxt) = 0;
        clusideas.cenfreq(z) = 0;
        clusideas.centime(z) = 0;
    end;
    if (y(maxf)-y(minf)) < 20
        tempmask(minf:maxf, mint:maxt) = 0;
        clusideas.cenfreq(z) = 0;
        clusideas.centime(z) = 0;
    end;
end;
clusideas.cenfreq = unique(clusideas.cenfreq);
clusideas.centime = unique(clusideas.centime);
