clc;
clear;
close all

%% Initialization 
set1 = reshape(table2array(readtable('Dataset_1 .csv')),[],[3]);
set2 = reshape(table2array(readtable('Dataset_2 .csv')),[],[3]);

alpha = 0.1;
minPts = 5;
MT = 100 ; %Merge Threshold

%% Start SOStream porccess
fprintf("SOStream on data set 1:\n")
tic;
[M,LastM,Deleted,Merged] = SOStreamProccess(set1,alpha,minPts,MT);
toc
Draw(M,LastM,Deleted,Merged);

fprintf("\n*********\n\n")
fprintf("SOStream on data set 2:\n")

tic;
[M,LastM,Deleted,Merged] = SOStreamProccess(set2,alpha,minPts,MT);
toc
Draw(M,LastM,Deleted,Merged);


%% SOStream Main Proccess function
function [M,LastM,Deleted,Merged] = SOStreamProccess(DS,alpha,minPts,MT)
M ={};
LastWin = [];
LastM = [];
Deleted=0;
Merged=0;

for iter=1:length(DS)
    vt = DS(iter,:);
    
    if iter == 1;LastM = [];
    else;LastM = M{end};end
    
    win = minDist(vt,LastM);
    
    if length(LastM) > minPts
        [winN,win.r] = findNeighbors(win,LastM,minPts);
        
        if dst(newCluster(vt), win) <= win.r
            [winN,win]=updateCluster(win,vt,alpha,winN);
        else
            nc = newCluster(vt);
        end
        
    overlap = findOverlap(win,winN);
    
    if ~isempty(overlap)
        
        [mc, dc] = mergeCluster(win, overlap, MT);
        Deleted = Deleted + length(dc);
        Merged = Merged + length(mc);
        for i=dc;LastM = deleteLastM(LastM,i);end
        LastM = [LastM mc];
        
        
    end
    
    else
        nc = newCluster(vt);
        LastM = [LastM nc];
    end
    M = [M {copyLastM(LastM)}];
    LastWin = win;
end
end

%% Draw Plots funciton
function Draw(M,LastM,Deleted,Merged)
fprintf("Total number of final clusters : %d\n",length(LastM));
fprintf("Total number of deleted clusters : %d\n",Deleted);
fprintf("Total number of merged clusters : %d\n",Merged);

dots = zeros(1,length(M));
figure
for i=1:length(M)
    dots(i) = length(M{1,i});
    subplot(2,1,1);
    
    t= title("Changes in the number of clusters");
    t.FontSize = 16;

    
    grid on;
    plot(i,length(M{1,i}),'r.')
    hold on 
    pause(0.00000001)
end
subplot(2,1,2);
plot([1:length(dots)],dots,'-x')
end

%% Other functions
function [mc, dc] = mergeCluster(win, ovc, MT)
    dc = [];
    mc = [];
    
    for Ni=ovc
        if dst(win,Ni)< MT
            if isempty(dc)
                dc = [dc win];
                mc = MicroCluster(win.n,win.r,win.c);
            end
            mc = merge(mc,Ni);
            dc = [dc Ni];
        end
    end
    
end

function nc = merge(ca, cb)
    %new cluster centroid
    ncc = weighted_mean(ca.c,cb.c,ca.n,cb.n);
    %new cluster radius
    ncr = dst(ca,cb) + max(ca.r,cb.r);
    %new cluster
    nc = MicroCluster((ca.n + cb.n), ncr, ncc);
end

function ov = findOverlap(win, winN)
    overlap = [];
    for N = winN
        if ~isequal(N,win)
            if (dst(win,N) - (win.r + N.r)) < 0
                overlap = [overlap N];
            end
        end
    end
    ov = overlap;
end

function out = weighted_mean(a,b,wa,wb)
    out = (a*wa + b*wb)/(wa + wb);
end
 
function [UwinN,Uwin]=updateCluster(win, vt, alpha,winN)
    win.c = weighted_mean(vt,win.c,1,win.n);
    win.n = win.n + 1;
    widthN = win.r^2;
    for cluster=winN
         influence = exp(-dst(cluster,win)/(2 * widthN));
         cluster.c=cluster.c+alpha*influence*(win.r-cluster.r);
    end
    UwinN = winN;
    Uwin = win;
end

function [nrs,r] = findNeighbors(win,M,minPts)
    if length(M) > minPts
        winDistN = zeros(1,length(M));
        for N=1:length(M)
            winDistN(N) = dst(win,M(N));
        end
        winDistN = sort(winDistN);
        %kdist: represent the radius (threshold) of the winning cluster
        kDist = winDistN(minPts);
        win.r = kDist;
        
        winNN =[];
        for i=1:length(winDistN)
            if winDistN(i) <= kDist
                winNN =[winNN M(i)];
            end
        end
        nrs = winNN;
        r= kDist;
    else
        nrs = [];
        r =0;
    end
end

function out = newCluster(c)
      out = MicroCluster(1,0,c);
%     out = {n,r,c}; % [number of point, radius ,centoid]
end

function d = dst(c1,c2)
    a = c1.c;
    b = c2.c;
    d = sqrt((a(1)-b(1))^2 + (a(2)-b(2))^2 + (a(3)-b(3))^2);
end

function m = minDist(vt,M)
    if ~isempty(M)
        mindst = inf;
        min=[];
        vtc = newCluster(vt);
        for i=1:size(M,2)
            if mindst > dst(vtc,M(i))
                min = M(i);
                mindst = dst(vtc,M(i));
            end
        end
        m = min;
    else
        m = [];
    end
end

function out = copyLastM(Last)
    out = [];
    for m=Last
        out = [out copy(m)];
    end
end

function Last = deleteLastM(Last,del)
    newLast = [];
    for c=Last
        if ~isequal(del,c)
            newLast = [newLast c];
        end
    end
    Last = newLast;
end

