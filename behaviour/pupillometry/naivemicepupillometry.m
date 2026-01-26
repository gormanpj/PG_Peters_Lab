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
dataDir = 'C:\Users\pgorman\Documents\SLEAP\Projects\thirdnaivepassivetest2';  
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

% Grab the third recording day and time for each mouse 

thirdPassives = cell(numel(mice),2); 

for m = 1:numel(mice) 
    passives = plab.find_recordings(mice(m), [], 'lcr_passive'); 
    thirdPassives(m,1) = cellstr(passives(3).day); 
    thirdPassives(m,2) = passives(3).recording(1); 
end 

frameStims = cell(numel(mice),3); 
for n = 1:numel(mice) 
    animal = mice{n}; 
    rec_day = thirdPassives{n,1}; 
    rec_time = thirdPassives{n,2}; 
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

diameterZAllFlip = cellfun(@transpose, diameterZ_all, 'UniformOutput', false);

preSec = 10.0;   % seconds before stimulus
postSec = 3.0;  % seconds after stimulus
orientations = [-90, 0, 90];    % expected stimulus orientations
nOr = numel(orientations);

% containers: for each orientation, a cell of trial matrices (each row = trial, cols = timepoints)
trials_by_ori = cell(nOr,1);

% loop videos
for vid = 1:numel(diameterZAllFlip)
    diam = diameterZAllFlip{vid};      % vector [frames x 1] in Z
    if isempty(diam), continue; end

    % get frame times and stim times/orientations for this video
    frameTimes = frameStims{vid,1};  % expected 1 x nFrames or vector
    stimTimes  = frameStims{vid,2};  % expected vector of stim timestamps
    stimOris   = frameStims{vid,3};  % expected vector of orientations

    % normalize shapes to row vectors
    frameTimes = frameTimes(:)';     % 1 x nFrameTimes
    stimTimes  = stimTimes(:)';      % 1 x nStimTimes
    stimOris   = stimOris(:)';       % 1 x nStimOris

    % if frameTimes longer than diam, drop tail frames
    if length(frameTimes) > length(diam)
        frameTimes = frameTimes(1:length(diam));
    end
    nFrames = length(frameTimes);

    % find frame sampling interval (median to be robust)
    if nFrames < 2
        warning('Video %d: too few frames', vid); continue;
    end
    frameDt = median(diff(frameTimes)); % seconds per frame

    % define peri-stim window in frames
    winPre  = round(preSec / frameDt);
    winPost = round(postSec / frameDt);
    winLen  = winPre + winPost + 1;
    tRel = (-winPre:winPost) * frameDt; % time vector in seconds for plotting

    % map each stimulus to nearest frame index (nearest frameTime)
    % use interp1 on frameTimes->indices for robustness
    frameIdxForStim = round(interp1(frameTimes, 1:nFrames, stimTimes, 'nearest', 'extrap'));

    nUse = length(stimTimes);
    % for each stimulus, extract peri-stim diameter segment and store by orientation
    for s = 1:nUse
        ori = stimOris(s);
        % skip orientations not in our list
        idxOri = find(orientations == ori, 1);
        if isempty(idxOri), continue; end

        centerIdx = frameIdxForStim(s);
        idxRange = (centerIdx - winPre) : (centerIdx + winPost);

        % create trial vector, pad out-of-bounds with NaN
        trial = nan(1, winLen);
        validMask = idxRange >= 1 & idxRange <= nFrames;
        if any(validMask)
            trial(validMask) = diam(idxRange(validMask));
        end

        % append to trials_by_ori{idxOri} as row
        if isempty(trials_by_ori{idxOri})
            trials_by_ori{idxOri} = trial;
        else
            trials_by_ori{idxOri}(end+1, :) = trial;
        end
    end
end

% Concatenate all trials across videos already occurred; now compute mean and SEM per orientation
psth_mean = cell(nOr,1);
psth_sem  = cell(nOr,1);
numTrials = zeros(nOr,1);

for k = 1:nOr
    T = trials_by_ori{k}; % trials x timepoints
    if isempty(T)
        psth_mean{k} = nan(1, winLen);
        psth_sem{k}  = nan(1, winLen);
        numTrials(k) = 0;
        continue;
    end
    numTrials(k) = size(T,1);
    psth_mean{k} = nanmean(T,1);                      % <--- mean computed here
    denom = sqrt(sum(~isnan(T),1));
    denom(denom==0) = NaN;
    psth_sem{k} = nanstd(T,0,1) ./ denom;             % <--- sem computed here
end

plotColors = [ ...
    0.00, 0.45, 0.74;   % -90 deg (blue)
    0.85, 0.33, 0.10;   %   0 deg (orange)
    0.47, 0.67, 0.19];  %  90 deg (green)

% --- Z-score sanity checks per video (prints mean & std of raw diameters used to compute Z) ---
fprintf('\nZ-score sanity check per video (mean, std of raw diameters used to compute Z):\n');
for vid = 1:numel(diameterPx_all)
    px = diameterPx_all{vid};
    if isempty(px), continue; end
    mu_px = nanmean(px);
    sigma_px = nanstd(px);
    fprintf('Video %d: mean=%.3f px, std=%.3f px\n', vid, mu_px, sigma_px);
end

% --- Plot PSTHs with shaded SEM ---
figure('Name','PSTHs by orientation (shaded SEM)', 'Color', 'w');
hold on;
hLine = gobjects(nOr,1);
hFill = gobjects(nOr,1);

for k = 1:nOr
    col = plotColors(k,:);
    mu = psth_mean{k}(:)';   % ensure row
    se = psth_sem{k}(:)';    % ensure row
    if all(isnan(mu))
        continue;
    end

    x = tRel(:)';            % row vector
    y_mean = mu;
    y_upper = mu + se;
    y_lower = mu - se;

    % Build polygon for shaded area
    xv = [x, fliplr(x)];
    yv = [y_upper, fliplr(y_lower)];

    % Draw shaded region (behind the mean line)
    hFill(k) = fill(xv, yv, col, 'FaceAlpha', 0.18, 'EdgeColor', 'none');

    % Draw mean line on top
    hLine(k) = plot(x, y_mean, 'LineWidth', 1.8, 'Color', col);
end

% cosmetics
xlabel('Seconds from stimulus onset');
ylabel('Pupil diameter (Z)');
legendTextAll = arrayfun(@(o,n) sprintf('%d° (n=%d)', o, n), orientations(:), numTrials(:), 'UniformOutput', false);

% Which orientations actually have PSTH data?
hasData = ~cellfun(@(m) all(isnan(m)), psth_mean);   % logical nOr x 1

% Filter handles and labels to only those plotted
validHandles = hLine(hasData);          % only plotted line handles
validLabels  = legendTextAll(hasData);  % corresponding labels

% Now create legend (only if at least one plotted line)
if any(hasData)
    legend(validHandles, validLabels, 'Location', 'best');
end
title(sprintf('PSTHs (pre=%.2fs post=%.2fs)', preSec, postSec));
grid on;
set(gca,'Box','off');

% --- Baseline checks: compute mean baseline (pre-stim window) for each orientation ---
fprintf('\nPSTH baseline means (should be ~0 if z-scoring OK):\n');
baselineIdx = find(tRel < 0); % indices for pre-stimulus times
for k = 1:nOr
    mu = psth_mean{k};
    if all(isnan(mu))
        fprintf('Ori %d: no data\n', orientations(k));
        continue;
    end
    baselineMean = nanmean(mu(baselineIdx));
    baselineSEM  = nanstd(mu(baselineIdx)) / sqrt(sum(~isnan(mu(baselineIdx))));
    fprintf('Ori %d: baseline mean = %.4f (SEM=%.4f), nTrials=%d\n', orientations(k), baselineMean, baselineSEM, numTrials(k));
end

hold off;