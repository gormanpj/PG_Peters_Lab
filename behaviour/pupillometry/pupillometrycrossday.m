% Script to plot pupil responses to lcr_passive stimuli across days, for
% one mouse

animalID = 'AM022';

% Use a folder
dataDir = 'C:\Users\pgorman\Documents\SLEAP\Outputs_Sorted\lcr_passive\AM022';  
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


%% Get recording times 

%% Aligning pupillometry to stimuli
% Follows on from processed pupil videos using hdf5 files (in same order as
% below)

% Grab the recording days and time for the mouse 

allPassives = cell(numel(fileList),2); 

passiveM = plab.find_recordings(animalID, [], 'lcr_passive')

% Now get frame times for each recording, 

[s,recs] = size(passiveM);

frameStims = cell(recs,3); 
for n = 1:recs 
    animal = passiveM(n).animal; 
    rec_day = passiveM(n).day; 
    rec_time = passiveM(n).recording{1}; 
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
diameterZAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterZAllFlip, 'UniformOutput', false);
diameterPxAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterPxAllFlip, 'UniformOutput', false);

% Get derivatives
diameterZAllFlipFiltDeriv = cellfun(@diff, diameterZAllFlipFilt, 'UniformOutput', false);
diameterPxAllFlipFiltDeriv = cellfun(@diff, diameterPxAllFlipFilt, 'UniformOutput', false);

%% Dilation analysis

% PSTH window
pre = 20;
post = 60;
x = -pre:post;                             % x-axis window
ncols = 4;                                 % columns in tiled layout
nrows = ceil(nFiles / ncols);

figure('Name','Dilation PSTHs','NumberTitle','off');
tiledlayout(nrows, ncols, 'TileSpacing','compact', 'Padding','compact');

for n = 1:nFiles
    nexttile;
    stimPositionsVec = frameStims{n,3};    % stimulus position vector 
    stimTimes = frameStims{n,2};           % stimulus onset times 
    frameTimes = frameStims{n,1};          % frame time vector 

    % Find stimulus onset times for each position 
    leftPres   = stimTimes(stimPositionsVec == -90);
    centerPres = stimTimes(stimPositionsVec == 0);
    rightPres  = stimTimes(stimPositionsVec == 90);

    % Convert stim times to frame indices 
    leftFrameIdx   = interp1(frameTimes, 1:numel(frameTimes), leftPres,   'previous', NaN)';
    centerFrameIdx = interp1(frameTimes, 1:numel(frameTimes), centerPres, 'previous', NaN)';
    rightFrameIdx  = interp1(frameTimes, 1:numel(frameTimes), rightPres,  'previous', NaN)';

    % Build trial x time matrices 
    leftRange   = leftFrameIdx   + (-pre:post);  
    centerRange = centerFrameIdx + (-pre:post);
    rightRange  = rightFrameIdx  + (-pre:post);

    leftData   = diameterZAllFlipFiltDeriv{n}(leftRange);    
    centerData = diameterZAllFlipFiltDeriv{n}(centerRange);
    rightData  = diameterZAllFlipFiltDeriv{n}(rightRange);

    % Compute mean and SEM across trials 
    meanLeft   = nanmean(leftData, 1);
    semLeft    = nanstd(leftData, 0, 1) ./ sqrt(sum(~isnan(leftData),1));

    meanCenter = nanmean(centerData, 1);
    semCenter  = nanstd(centerData, 0, 1) ./ sqrt(sum(~isnan(centerData),1));

    meanRight  = nanmean(rightData, 1);
    semRight   = nanstd(rightData, 0, 1) ./ sqrt(sum(~isnan(rightData),1));
    recAnimal = passiveM(n).animal;   % the animal ID for this tile
    recDay    = passiveM(n).day;      % recording day string, e.g. '2023-12-08'
    
    % Use aniDates from before to label learning day
    matchIdx = find(strcmp(aniDates(:,1), recAnimal), 1);
    learningLabel = '';
    if ~isempty(matchIdx)
        learnDateVal = aniDates{matchIdx, 2};
        % treat numeric NaN as missing; if a string present, compare
        if ischar(learnDateVal) || isstring(learnDateVal)
            if strcmp(string(learnDateVal), string(recDay))
                learningLabel = ' (learning day)';
            end
        end
    end
    
    % set title
    recAnimal = passiveM(n).animal;   
    recDay    = passiveM(n).day;      
    title(sprintf('%s — %s%s', animalID, recDay, learningLabel), ...
          'Interpreter','none', 'FontSize', 9);

    % Plot mean + SEM fills
    hold on;
    % Left (blue)
    xl = [x, fliplr(x)];
    yl = [meanLeft + semLeft, fliplr(meanLeft - semLeft)];
    fill(xl, yl, 'b', 'EdgeColor','none', 'FaceAlpha', 0.2);
    plot(x, meanLeft, 'b', 'LineWidth', 1.2);

    % Center (black)
    xc = xl;
    yc = [meanCenter + semCenter, fliplr(meanCenter - semCenter)];
    fill(xc, yc, 'k', 'EdgeColor','none', 'FaceAlpha', 0.2);
    plot(x, meanCenter, 'k', 'LineWidth', 1.2);

    % Right (red)
    xr = xl;
    yr = [meanRight + semRight, fliplr(meanRight - semRight)];
    fill(xr, yr, 'r', 'EdgeColor','none', 'FaceAlpha', 0.2);
    plot(x, meanRight, 'r', 'LineWidth', 1.2);

    xline(0,'--');
    yline(0,'Color','k','Alpha',0.5);
    xlabel('Frames (relative to stimulus)');
    ylabel('Pupil diameter derivative (arb)');
    ylim([-0.05, 0.05]);
   
end
