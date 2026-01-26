% Proof of principle for pupil diameter tracking (getting labels straight from
% HDF5 files)

% Next steps will be to align to stimuli (load data, maybe this script) and
% then integrate the model and prediction process into matlab itself


% Manually specify files
% fileList = {
% %    'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.000_mousecam.analysis.h5';
% %    'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.001_mousecam.analysis.h5';
%     'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.002_mousecam.analysis.h5';
%     'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.003_mousecam.analysis.h5';
%     'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.004_mousecam.analysis.h5';
%     'C:\Users\pgorman\Documents\SLEAP\Projects\testsetAM011\testsetAM011.005_mousecam.analysis.h5'
% };

% Alt: use a folder
dataDir = 'C:\Users\pgorman\Documents\SLEAP\Projects\AMfinaltrainedpassivetest2';  
pattern = fullfile(dataDir, '*.analysis.h5');  
files = dir(pattern);                          
fileList = fullfile({files.folder}, {files.name})';

scoreThresh = 0.6;    % confidence threshold for SLEAP point scores
nFiles = numel(fileList);

% Preallocate outputs. Cells stop lengths from being an issue.
diameterPx_all   = cell(nFiles,1);
diameterZ_all    = cell(nFiles,1);
radius_all       = cell(nFiles,1);
center_all       = cell(nFiles,1);
fitRmse_all      = cell(nFiles,1);

% Big loop to get the above for each video

for idx = 1:nFiles

    h5file = fileList{idx};

    % Read datasets

    t = h5read(h5file,'/tracks');       % frames x nodes x 2
    pointScores = h5read(h5file,'/point_scores')'; % frames x nodes
    instanceScores = h5read(h5file, '/instance_scores')'; % transpose to 1 x frames

    sz = size(t);
    numFrames = sz(1);
    numNodes= sz(2);

    X = squeeze(t(:,:,1))'; % becomes [numFrames x numNodes]
    Y = squeeze(t(:,:,2))'; % [numFrames x numNodes]

    % Mask low-confidence instances 
    X(:,  instanceScores< scoreThresh) = NaN;
    Y(:,  instanceScores< scoreThresh) = NaN;
    
    fprintf('%d bad instances dropped, %g3%% of video \n', sum(sum(instanceScores<scoreThresh)), ((sum(sum(instanceScores<scoreThresh)))./numFrames)*100);

    % Mask low-confidence points
    X(pointScores < scoreThresh) = NaN;
    Y(pointScores < scoreThresh) = NaN;

    fprintf('%d additional bad nodes dropped \n \n', sum(sum(pointScores<scoreThresh)));

    % Commented out: code to interpolate missing points. this is more
    % likely to cause major errors in tracking than to help. 
    % for n=1:numNodes
    %     xn = X(n,:);
    %     yn = Y(n,:);
    % 
    %     % linear interpolation for interior NaNs
    %     xn = fillmissing(xn,'linear','EndValues','nearest');
    %     yn = fillmissing(yn,'linear','EndValues','nearest');
    % 
    %     X(n,:) = xn;
    %     Y(n,:) = yn;
    % end

    % Fit circles to frames
    radius = nan(numFrames,1);
    center = nan(numFrames,2);
    diameterPx = nan(numFrames,1);
    fitRmse = nan(numFrames,1);

    for f = 1:numFrames
        xpts = X(:,f);
        ypts = Y(:,f);
        valid = ~isnan(xpts) & ~isnan(ypts);
        if nnz(valid) < 3
            % Not enough points to fit a circle
            radius(f) = NaN;
            center(f,:) = [NaN NaN];
            diameterPx(f) = NaN;
            fitRmse(f) = NaN;
            continue;
        end

        xg = xpts(valid);
        yg = ypts(valid);
    
        % Algebraic least-squares fit:
        % Solve for a,b,c in x^2 + y^2 + a*x + b*y + c = 0
        A = [xg, yg, ones(length(xg),1)];
        bvec = -(xg.^2 + yg.^2);
    
        % Solve linear system (least squares)
        p = A \ bvec;        % p = [a; b; c]
        a = p(1); bpar = p(2); c = p(3);
    
        xc = -a/2;
        yc = -bpar/2;
        radTerm = (a^2 + bpar^2)/4 - c;
        if radTerm <= 0
            % numerical degeneracy -> treat as invalid
            radius(f) = NaN;
            center(f,:) = [NaN NaN];
            diameterPx(f) = NaN;
            fitRmse(f) = NaN;
            continue;
        end
        R = sqrt(radTerm);
    
        % Store results
        radius(f) = R;
        center(f,:) = [xc yc];
        diameterPx(f) = 2*R;
    
        % Compute RMSE of radial residuals as a fit quality metric. Might have
        % to fit ovals instead if the fit is bad
        dists = hypot(xg - xc, yg - yc);
        residuals = dists - R;
        fitRmse(f) = sqrt(mean(residuals.^2));
    end

    % Now Z-score
    mu = nanmean(diameterPx);
    sig = nanstd(diameterPx);
    if sig == 0 || isnan(sig)
        diameterZ = nan(size(diameterPx));
    else
        diameterZ = (diameterPx - mu) ./ sig;
    end
    
    % Store results
    diameterPx_all{idx} = diameterPx;
    diameterZ_all{idx}  = diameterZ;
    radius_all{idx}     = radius;
    center_all{idx}     = center;
    fitRmse_all{idx}    = fitRmse;

end

% Plot for raw diameters
figure('Name','Raw pupil diameters (px, overlay)');
hold on;
colors = lines(nFiles);
for idx = 1:nFiles
    d = diameterPx_all{idx};
    if isempty(d), continue; end
    tvec = 1:length(d);
    plot(tvec, d, 'Color', colors(mod(idx-1,size(colors,1))+1,:), 'DisplayName', sprintf('file %d', idx));
end
xlabel('Frame'); ylabel('Pupil diameter (px)');
legend('show');
title('Raw pupil diameters (per-file overlay)');
hold off;

% Plot for z-scored diameters
figure('Name','Z-scored pupil diameters (per-file overlay)');
hold on;
colors = lines(nFiles);
for idx = 1:nFiles
    z = diameterZ_all{idx};
    if isempty(z), continue; end
    tvec = 1:length(z);
    plot(tvec, z, 'Color', colors(mod(idx-1,size(colors,1))+1,:), 'DisplayName', sprintf('file %d', idx));
end
xlabel('Frame'); ylabel('Pupil diameter (z-score)');
legend('show');
title('Per-video z-scored pupil diameter (no temporal alignment)');
hold off;

% Print the RMSE values in case circular fitting isn't working well
fprintf('\nMean RMSE per file (px):\n');
for idx = 1:nFiles
    rm = fitRmse_all{idx};
    if isempty(rm) || all(isnan(rm))
        fprintf('File %d: no valid RMSE data\n', idx);
    else
        meanRmse = mean(rm(~isnan(rm)));
        fprintf('File %d: %.3f px\n', idx, meanRmse);
    end
end

%% Aligning pupillometry to stimuli
% Follows on from processed pupil videos using hdf5 files (in same order as
% below)

% All mice in naive group 
mice = {'AM011'; 
        'AM012'; 
        'AM014'; 
        'AM015'; 
        'AM016'; 
        'AM017'; 
        'AM018'; 
        'AM019'; 
        'AM021'; 
        'AM022'; 
        'AM026'; 
        'AM029'; 
        'AP023'; 
        'AP025'; 
        }; 

% Grab the final recording day and time for each mouse 

finalPassives = cell(numel(mice),2); 

for m = 1:numel(mice) 
    passives = plab.find_recordings(mice(m), [], 'lcr_passive'); 
    finalPassives(m,1) = cellstr(passives(end).day); 
    finalPassives(m,2) = passives(end).recording(end); 
end 

frameStims = cell(numel(mice),3); 
for n = 1:numel(mice) 
    animal = mice{n}; 
    rec_day = finalPassives{n,1}; 
    rec_time = finalPassives{n,2}; 
    verbose = false; % this turns on/off progress display in command line 

    % Loading components separately. Could maybe be more efficient but is
    % better than calling ap.load_recording
    ap.load_timelite;
    ap.load_mousecam;
    ap.load_bonsai; 
    
    % Grab the frame times and put them in one cell 
    frameStims(n,1) = mat2cell(mousecam_times,1); 
    frameStims(n,2) = mat2cell(stimOn_times',1); 
    frameStims(n,3) = mat2cell((vertcat(trial_events.values.TrialStimX))',1); 
end 

diameterPxAllFlip = cellfun(@transpose, diameterPx_all, 'UniformOutput', false);
diameterZAllFlip = cellfun(@transpose, diameterZ_all, 'UniformOutput', false);


% Now we have frameStims with mousecam times, stimOn_times, and stim
% orientations. we have diameterZAllFlip with z-scored pupil diameters in
% the same orientation

% Fourier analysis of pupil diameter signal

figure; tiledlayout(4,4); 
for r = 1:nFiles
    nexttile;
    ex = diameterPx_all{r,1};
    exfill = fillmissing(ex,'movmean', 120);
    % figure; hold on; % Allows us to glance for NaNs and weirdness
    % plot(exfill, 'Color', 'r'); 
    % plot(ex, 'Color', 'k'); 
    % hold off;
    
    fex = fft(exfill - mean(exfill)); % fast fourier transform
    T = 1/30;
    fq = 1/T;
    fourex = (0:length(exfill)-1)*fq/length(exfill);
    
    % plot(fourex, abs(fex));
    % xlabel('Freq (hz)'); ylabel('Magnitude');
    % title('Fourier transform of pupil diameter');
    
    n = length(exfill);
    fshift = (-n/2:n/2-1)*(fq/n);
    fexshift = fftshift(fex);
    plot(fshift,abs(fexshift))
    xlim([-15,-2]);
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title("video "+r)
end
% Getting 0-values (DC signal?) but also something like a peak at ~9.2hz.
% Weird

% Lowpass filter to get rid of jitteriness
diameterZAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4.5, 30), diameterZAllFlip, 'UniformOutput', false);
diameterPxAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4.5, 30), diameterPxAllFlip, 'UniformOutput', false);

% To plot this as derivatives

diameterZAllFlipFiltDeriv = cellfun(@diff, diameterZAllFlipFilt, 'UniformOutput', false);
diameterPxAllFlipFiltDeriv = cellfun(@diff, diameterPxAllFlipFilt, 'UniformOutput', false);

% Excluding movement 

% Movement artifacts
% Does the pupil dilation correspond to general activity of wheel movement?

% We'll have to take timelite.timestamps and wheel_move for each video to
% create a mask and use the mask to NaN out diameter traces
movMask = cell(nFiles, 1);

for m = 1:nFiles
    animal = mice{m}; 
    rec_day = finalPassives{m,1}; 
    rec_time = finalPassives{m,2}; 
    verbose = false; % this turns on/off progress display in command line 

    % Loading components separately. Could maybe be more efficient but is
    % better than calling ap.load_recording
    ap.load_timelite;
    ap.load_mousecam;

    ft = mousecam_times(:);
    t  = timelite.timestamps(:);
    move = wheel_move(:);
    
    edges = [ft(1) - (ft(2)-ft(1))/2 ; (ft(1:end-1)+ft(2:end))/2 ; ft(end) + (ft(end)-ft(end-1))/2];
    
    % Bin high-rate timestamps to frame intervals
    bin = discretize(t, edges);
    
    % Vectorized per-frame OR ("any movement in frame")
    valid = ~isnan(bin);  % only timestamps within the frame range
    movement_per_frame = accumarray(bin(valid), move(valid), [numel(ft), 1], @max, false);
    
    movMask{m} = movement_per_frame;
end

% Perfect fixed code that no longer kills matlab 
movDiameterPxAllFilt = diameterPxAllFlipFiltDeriv;

for n = 1:nFiles
    mask = movMask{n};                 % logical mask the same size as the data
    movDiameterPxAllFilt{n}(mask) = NaN;
end

%% Dilation analysis

%First let's just do a PSTH of pupil dilation around the stimuli
pre = 20;
post = 60;

leftMeans = [];
centerMeans = [];
rightMeans = [];

for n = 1:nFiles
    % Getting right, center and left frames from demo script and frameStims

    stimPositionsVec = frameStims{n,3}; % Grab the cells of all stimulus positions
    
    leftPres = frameStims{n,2}(find(stimPositionsVec == -90)); 
    centerPres = frameStims{n,2}(find(stimPositionsVec == 0));
    rightPres = frameStims{n,2}(find(stimPositionsVec == 90)); % Apply this order to the timelite array
    
    leftFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), leftPres, 'previous', NaN)';
    centerFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), centerPres, 'previous', NaN)';
    rightFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), rightPres, 'previous', NaN)';
    
    leftRange = leftFrameIdx + (-pre:post); 
    centerRange = centerFrameIdx + (-pre:post);
    rightRange = rightFrameIdx + (-pre:post); 
    
    meanLeftAngleRange= nanmean(movDiameterPxAllFilt{n}(leftRange), 1);
    meanCenterAngleRange= nanmean(movDiameterPxAllFilt{n}(centerRange), 1);
    meanRightAngleRange= nanmean(movDiameterPxAllFilt{n}(rightRange), 1);

    leftMeans(n,:) = meanLeftAngleRange;
    centerMeans(n,:) = meanCenterAngleRange;
    rightMeans(n,:) = meanRightAngleRange;
end

meanL = nanmean(leftMeans, 1);
meanC = nanmean(centerMeans, 1);
meanR = nanmean(rightMeans, 1);

semL = nanstd(leftMeans, 0, 1) ./ sqrt(sum(~isnan(leftMeans), 1));
semC = nanstd(centerMeans, 0, 1) ./ sqrt(sum(~isnan(centerMeans), 1));
semR = nanstd(rightMeans, 0, 1) ./ sqrt(sum(~isnan(rightMeans), 1));

x = -pre:post;

figure; hold on;

xl = [x, fliplr(x)];
yl = [meanL + semL, fliplr(meanL - semL)];
hL = fill(xl, yl, 'b', 'EdgeColor','none', 'FaceAlpha', 0.2);
plot(x, meanL, 'Color','b'); 
% for l = 1:nFiles
%     plot(x, leftMeans(l,:), '--', 'Color','b');
% end

xc = xl;
yc = [meanC + semC, fliplr(meanC - semC)];
hC = fill(xc, yc, 'k', 'EdgeColor','none', 'FaceAlpha', 0.2);
plot(x, meanC, 'Color','k'); 
% for c = 1:nFiles
%     plot(x, centerMeans(c,:), '--', 'Color','k');
% end

xr = xl;
yr = [meanR + semR, fliplr(meanR - semR)];
hR = fill(xr, yr, 'r', 'EdgeColor','none', 'FaceAlpha', 0.2);
plot(x, meanR, 'Color','r'); 
% for r = 1:nFiles
%     plot(x, rightMeans(r,:), '--', 'Color','r');
% end

xline(0,'--')
yline(0,'Color', 'k', 'Alpha',0.5);
ylim([-0.025, 0.025])
xlabel("Frames"); ylabel("Pupil diameter deriv. (arb)");
legend({'Left SEM','Left Mean','Center SEM','Center mean','Right SEM','Right mean'}, 'Location','best');
hold off;


% Sanity checks for dilation analysis
% pick a signal to visualize 
figure; tiledlayout(4,4);
for t = 1:nFiles
    nexttile;
    plot(diameterZAllFlipFilt{t});
    title("video "+ t)
end

% make an average video. This is a bit janky but I'm grabbing it from my plab demo code
sampleVideo = 5;
spv = frameStims{sampleVideo,3};
rp = frameStims{sampleVideo,2}(find(spv == 90));
rfi = interp1(frameStims{sampleVideo,1}, 1:numel(frameStims{sampleVideo,1}), rp, 'previous', NaN)';
rr = rfi + (-20:30); 

animal = mice{sampleVideo}; 
rec_day = finalPassives{sampleVideo,1}; 
rec_time = finalPassives{sampleVideo,2};
ap.load_mousecam;

vid = VideoReader(mousecam_fn); 
[nPres, nWindow] = size(rr); % Need the numbers of presentations and window size
[dim1,dim2] = size(read(vid,1));
firstframe = [size(read(vid,1)), nWindow];

avgRight = zeros(firstframe, "double"); 
for i = 1:nWindow
    frameIds = [rr(:,i)]; % Get the frame indices for column in question
    sumFrames = zeros(dim1, dim2, 'double');
    for j = 1:numel(frameIds) % Stack em
        frameIdx = frameIds(j);
        frame = read(vid, frameIdx);
        sumFrames = sumFrames + double(frame);
    end
    avgRight(:,:,i) = sumFrames /numel(frameIds); % Average pixel intensities into the avgRight matrix
end

averageRight = zeros(firstframe, "double");
for i = 1:nPres
    divPres = double(squeeze(read(vid, rr(i, [1, end]))))/nPres;
    averageRight = averageRight + double(divPres);
    %disp(i);
end

ap.imscroll(averageRight);


% Movement artifacts
% Does the pupil dilation correspond to general activity of wheel movement?

% We'll have to take timelite.timestamps and wheel_move for each video to
% create a mask and use the mask to NaN out diameter traces
movMask = cell(nFiles, 1);

for m = 1:nFiles
    animal = mice{m}; 
    rec_day = finalPassives{m,1}; 
    rec_time = finalPassives{m,2}; 
    verbose = false; % this turns on/off progress display in command line 

    % Loading components separately. Could maybe be more efficient but is
    % better than calling ap.load_recording
    ap.load_timelite;
    ap.load_mousecam;

    ft = mousecam_times(:);
    t  = timelite.timestamps(:);
    move = wheel_move(:);
    
    edges = [ft(1) - (ft(2)-ft(1))/2 ; (ft(1:end-1)+ft(2:end))/2 ; ft(end) + (ft(end)-ft(end-1))/2];
    
    % Bin high-rate timestamps to frame intervals
    bin = discretize(t, edges);
    
    % Vectorized per-frame OR ("any movement in frame")
    valid = ~isnan(bin);  % only timestamps within the frame range
    movement_per_frame = accumarray(bin(valid), move(valid), [numel(ft), 1], @max, false);
    
    movMask{m} = movement_per_frame;
end


%% Saccade analysis

% Find midpoint for each video

center_allAvgs = cellfun(@nanmean, center_all, 'UniformOutput', false);

%first scatters for comparison and troubleshooting
figure; tiledlayout(4,4);
for f = 1:nFiles
   
    nexttile;
    hold on;
    scatter(center_all{f,1}(:,1), center_all{f,1}(:,2), 2.5);
    title("video "+f);
    plot(center_allAvgs{f}(1), center_allAvgs{f}(2), '.', 'MarkerSize', 15, 'Color', 'r');
    hold off;
end

center_allAvgs{1}(2)

% as heatmaps
nBins = 50;
figure; tiledlayout(4,4);
for f = 1:nFiles
    nexttile
    [ct,Xedg,Yedg] = histcounts2(center_all{f,1}(:,2), center_all{f,1}(:,1),  nBins); 
    imagesc(Yedg,Xedg,ct);
    set(gca,'YDir','normal'); % why does imagesc do this
    title("video "+f)
end

% Now some sort of measure of deflection 
% I think we'll want to store vectors of direction and displacement for
% each frame

displacementAll = cell(nFiles, 1);
deflectionAll = cell(nFiles, 1);
angleDegAll = cell(nFiles, 1);
angleRadAll = cell(nFiles, 1);

for v = 1:nFiles
    vid = center_all{v,1};
    vidcenter = center_allAvgs{v};
    x = vidcenter(1);
    y = vidcenter(2);

    % displacement
    dx = vid(:,1) - x;
    dy = vid(:,2) - y;

    %distance
    deflection = hypot(dx,dy);
    angledeg = atan2d(dy,dx); % y first. learned the hard way
    anglerad = atan2(dy,dx); 

    %store
    displacementAll{v} = [dx,dy];
    deflectionAll{v} = deflection;
    angleDegAll{v} = angledeg;
    angleRadAll{v} = anglerad;
end

% Now to plot this. Let's just figure out one video at a time and then we
% can take averages

% Plot raw deflection
figure;
tiledlayout(6,5,'TileIndexing', 'columnmajor')
for d = 1:nFiles
    nexttile;
    plot(displacementAll{d}(:,1));
    ylabel("deflection(px)");
    xlabel("frames");
    ylim([-5,5]);
    title("X video "+d);
    nexttile;
    plot(displacementAll{d}(:,2), 'Color','k');
    ylabel("deflection(px)");
    xlabel("frames");
    ylim([-5,5]);
    title("Y video "+d);
end

figure; tiledlayout(4,4);
for f = 1:nFiles
    nexttile;
    hold on;
    scatter(displacementAll{f,1}(:,1), displacementAll{f,1}(:,2), 2.5);
    title("video "+f);
    plot(0, 0, '.', 'MarkerSize', 15, 'Color', 'r');
    hold off;
end

% Plot position (degrees)

figure;
tiledlayout(2,1,'TileIndexing', 'columnmajor')
for d = 1:1
    nexttile;
    plot(angleDegAll{d}(:,1));
    ylabel("position (degree)");
    xlabel("frames");
    title("Angle video "+d);
end

figure; tiledlayout(4,4)
for a = 1:nFiles
    nexttile
    polarhistogram(angleRadAll{a});
    title("video "+a)
end

% What's the best plot here? Like a heatmap of all vectors (plot displacement
% as points in bins) relative to center for all videos? Actually yes that's great. 

% We could normalize using average pupil diameter
%avgDiameterPxAll = cellfun(@nanmean, diameterPx_all);

% This normalizes using the longest distance
displacementNorm = cell(nFiles,1);
for d =1:nFiles
    displacementNorm{d,1} = displacementAll{d,1}/max(deflectionAll{d});
end

% Plot non-normalized vs normalized 
figure; hold on;
for s = 1:nFiles
    scatter(displacementAll{s,1}(:,1), displacementAll{s,1}(:,2), 2.5);
end
hold off;

figure; hold on;
for s = 1:nFiles
    scatter(displacementNorm{s,1}(:,1), displacementNorm{s,1}(:,2), 2.5);
end
hold off;

% General plot - something isn't right here

nBins = 100;
Xedg = linspace(-1,1,nBins+1);
Yedg = linspace(-1,1,nBins+1);
ctTotal = zeros(nBins);
xCenters = Xedg(1:end-1) + diff(Xedg)/2;
yCenters = Yedg(1:end-1) + diff(Yedg)/2;

for h = 1:nFiles
    [ct,Xedg,Yedg] = histcounts2(displacementNorm{h}(:,2), displacementNorm{h}(:,1),  Xedg,Yedg);
    ctTotal = ctTotal + ct;
end

figure; hold on;
imagesc(yCenters,xCenters,ctTotal);
%set(gca,'YDir','normal'); % why does imagesc do this
plot(0, 0, '.', 'MarkerSize', 15, 'Color', 'r');
colormap parula;
title('Normalized displacement histogram, all videos')
axis tight;
hold off;

%figure; scatter(ctTotal, 2.5);

% How does this change for stimulus presentations?

% lets make something like a measure of probability relative to total trial
alpha = 1; % tune: 0.1..5 depending on counts
pTotal = (ctTotal + alpha) / sum(ctTotal(:) + alpha);

%nBins = 16;
nFramesAfterStim = 40;        
stimOri = [-90, 0, 90];      
nOri = numel(stimOri);

ctOri = cell(1, nOri);
for o=1:nOri
    ctOri{o} = zeros(nBins, nBins);
end

for vid = 1:nFiles
    frameTimestamps = frameStims{vid, 1};   
    stimTimes       = frameStims{vid, 2};    
    stimOris        = frameStims{vid, 3};    
    
    d = displacementNorm{vid};   % Nx2 (X, Y)
 
    nFrames = size(d,1);
    
    % For each orientation, find stimulus presentations and gather frame indices
    for o = 1:nOri
        oriVal = stimOri(o);
        stimIdxs = find(stimOris == oriVal);   
        
        stimTimesThisOri = stimTimes(stimIdxs);
        
        frameIdx = interp1(frameTimestamps, 1:numel(frameTimestamps), stimTimesThisOri, 'previous', NaN)';
        
        % Collect displacement samples for the 40 frames after each stimulus
        collected = []; % will store [X Y] rows
        for k = 1:numel(frameIdx)
            fidx = frameIdx(k);
            if isnan(fidx)
                continue
            end
            % window: next frame after stimulus up to nFramesAfterStim frames after
            winStart = fidx + 1;
            winEnd   = fidx + nFramesAfterStim;
            % clamp to available frames
            winStart = max(1, winStart);
            winEnd   = min(nFrames, winEnd);
            if winEnd >= winStart
                collected = [collected; d(winStart:winEnd, :)]; %#ok<AGROW>
            end
        end
        
        
        [ct, Xedg, Yedg] = histcounts2(collected(:,2), collected(:,1), Xedg,Yedg);
        ctOri{o} = ctOri{o} + ct;
    end
end

% probs
pOri = cell(1,nOri);
for o = 1:nOri
    pOri{o} = (ctOri{o} + alpha) / sum(ctOri{o}(:) + alpha);
end

epsVal = 1e-12;
ratioOri = cell(1,nOri);
logRatioOri = cell(1,nOri);
for o = 1:nOri
    ratioOri{o} = (pOri{o}) ./ (pTotal + epsVal);
    logRatioOri{o} = log2(ratioOri{o} + epsVal);
end

% smoothing (optional)
smoothSigma = 0.8; % tune
smoothedLog = cell(1,nOri);
smoothedProb = cell(1,nOri);
for o = 1:nOri
    smoothedLog{o} = imgaussfilt(logRatioOri{o}, smoothSigma);
    smoothedProb{o} = imgaussfilt(pOri{o}, smoothSigma);
end

figure;
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

for o = 1:nOri
    % 1) raw probability map (smoothed)
    nexttile
    imagesc(Yedg, Xedg, smoothedProb{o});
    axis xy; set(gca,'YDir','normal'); hold on;
    plot(0,0,'.','MarkerSize',18,'Color','k');
    colorbar;
    title(sprintf('probability of displacement relative to whole trial, ori=%d°', stimOri(o)));
    xlabel('X displacement (norm)'); ylabel('Y displacement (norm)');
    axis tight; hold off;
end





% This is an attempt at graphing angle in response to the three stimulus
% orientations. Currently its uninformative. Also the angle averaging is
% wrong
leftMeans = [];
centerMeans = [];
rightMeans = [];

for n = 1:nFiles
    % Getting right, center and left frames from demo script and frameStims

    stimPositionsVec = frameStims{n,3}; % Grab the cells of all stimulus positions
    
    leftPres = frameStims{n,2}(find(stimPositionsVec == -90)); 
    centerPres = frameStims{n,2}(find(stimPositionsVec == 0));
    rightPres = frameStims{n,2}(find(stimPositionsVec == 90)); % Apply this order to the timelite array
    
    leftFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), leftPres, 'previous', NaN)';
    centerFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), centerPres, 'previous', NaN)';
    rightFrameIdx = interp1(frameStims{n,1}, 1:numel(frameStims{n,1}), rightPres, 'previous', NaN)';
    
    leftRange = leftFrameIdx + (-20:20); 
    centerRange = centerFrameIdx + (-20:20);
    rightRange = rightFrameIdx + (-20:20); 
    
    meanLeftAngleRange= mean(angleAll{n}(leftRange), 1);
    meanCenterAngleRange= mean(angleAll{n}(centerRange), 1);
    meanRightAngleRange= mean(angleAll{n}(rightRange), 1);

    leftMeans(n,:) = meanLeftAngleRange;
    centerMeans(n,:) = meanCenterAngleRange;
    rightMeans(n,:) = meanRightAngleRange;
end

figure; plot(nanmean(leftMeans))
figure; plot(nanmean(centerMeans))
figure; plot(nanmean(rightMeans))