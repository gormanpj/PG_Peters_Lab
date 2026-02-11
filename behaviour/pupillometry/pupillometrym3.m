%% Load behavior, add 'learned_days' and 'days_from_learning' fields
data_path = fullfile(plab.locations.server_path,'Lab','Papers','Marica_2025','data');
% Set stat and p-value to define learning day
use_stat = 'firstmove_mean';
learn_stat_p = 0.05;
% Load behavior
load(fullfile(data_path,'bhv'));
% Set "learned_days" and "days_from_learning"
bhv.learned_days = cellfun(@(x) x < learn_stat_p,bhv.(['stimwheel_pval_',use_stat]));
for curr_animal = unique(bhv.animal)'
    curr_rows = strcmp(bhv.animal,curr_animal);
    bhv.days_from_learning(curr_rows) = ...
        (1:sum(curr_rows))' - ...
        max([NaN,find(bhv.learned_days(curr_rows),1)]);
end

%% Grab the learning day for each animal
% See unique animals in bhv
aniUniq = unique(bhv.animal);
% initialize cell to store learning days
aniDates = cell(numel(aniUniq),2);
for n = 1:numel(aniUniq); 
    id = aniUniq{n};
    aniDates(n,1) = {id};
    % Get only recordings for this animal
    daysID = bhv(strcmp(bhv.animal, id),:);
    if any(daysID.days_from_learning == 0) % check if there was a learning day
        ld = daysID(daysID.days_from_learning == 0, :);
        aniDates(n,2) = ld.rec_day; % set learning date
    else
        aniDates(n,2) = {NaN}; % if no learning date, NaN
    end
end

%% Script to plot pupil responses to lcr_passive stimuli across days, for many mice

mouseIDs = unique(bhv.animal); % grabbing mice in question from bhv table

scoreThresh = 0.6;    % confidence threshold for SLEAP point scores

% PSTH window
pre = 20;
post = 60;
x = -pre:post; 

mouseLabels = mouseIDs;

% pooled outputs 
psthData_all       = [];
orientationIdx_all = [];
learnDayIdx_all    = [];
mouseIdx_all       = [];    

% Loop across all mice
for m = 1:numel(mouseIDs)
    animalID = mouseIDs{m};
    fprintf(animalID); fprintf(' \n \n');
    dataDir = fullfile('C:\Users\pgorman\Documents\SLEAP\Outputs_Sorted\lcr_passive', animalID);  
    pattern = fullfile(dataDir, '*.analysis.h5');  
    files = dir(pattern);                          
    fileList = fullfile({files.folder}, {files.name})';
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
        
        fprintf('%d bad instances dropped, %g3%% of video.', sum(sum(instanceScores<scoreThresh)), ((sum(sum(instanceScores<scoreThresh)))./numFrames)*100);
    
        % Mask low-confidence points
        X(pointScores < scoreThresh) = NaN;
        Y(pointScores < scoreThresh) = NaN;
    
        fprintf(' %d additional bad nodes dropped \n \n', sum(sum(pointScores<scoreThresh)));
    
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
            
            % Replace here with median difference b/t opposed points 

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

    %% Aligning pupillometry to stimuli
    % Follows on from processed pupil videos using hdf5 files (in same order as
    % below)
    
    % Grab the recording days and time for the mouse 
    
    allPassives = cell(numel(fileList),2); 
    
    passiveM = plab.find_recordings(animalID, [], 'lcr_passive');
    
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

    % Lowpass filter to get rid of jitteriness
    diameterZAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterZAllFlip, 'UniformOutput', false);
    diameterPxAllFlipFilt = cellfun(@(x) lowpass(fillmissing(x, "linear"), 4, 30), diameterPxAllFlip, 'UniformOutput', false);

    % try movemean after
    diameterZAllFlipFiltMov = cellfun(@(x) movmean(x, [6 0]), diameterZAllFlipFilt, 'UniformOutput', false);
    diameterPxAllFlipFiltMov = cellfun(@(x) movmean(x, [6 0]), diameterPxAllFlipFilt, 'UniformOutput', false);

    % Get derivatives
    diameterZAllFlipFiltDeriv = cellfun(@diff, diameterZAllFlipFiltMov, 'UniformOutput', false);
    diameterPxAllFlipFiltDeriv = cellfun(@diff, diameterPxAllFlipFiltMov, 'UniformOutput', false);

    % Output to use
    pupilPerFile = diameterZAllFlipFiltDeriv;  
    
    % Prepare accumulators
    psthData       = [];   % will become [nTrials x (pre+post+1)]
    orientationIdx = [];   % will become [nTrials x 1] (1=left,2=center,3=right)
    dayIdx         = [];   % will become [nTrials x 1] (recording index n)
    orientationLab = {};   
    recDayLabels   = {};   
    
    nWindow = numel(-pre:post);
    nRecs = min(nFiles, size(frameStims,1)); 
    
    for n = 1:nRecs
        % grab stimulus/frame info for this recording 
        frameTimes = frameStims{n,1};
        stimTimes  = frameStims{n,2};
        stimPosVec = frameStims{n,3};    
        
        % pupil trace for this recording
        trace = pupilPerFile{n}(:); % ensure column vector
        nFrames = numel(trace);
        
        % convert stim times to nearest (previous) frame indices
        % ensure stimTimes and stimPosVec are column vectors 
        stimTimes = frameStims{n,2};
        stimPosVec = frameStims{n,3};
        stimTimes = stimTimes(:);
        stimPosVec = stimPosVec(:);
    
        stimFrameIdx = interp1(frameTimes, (1:numel(frameTimes))', stimTimes, 'previous', NaN);
    
        % if lengths differ for some reason, trim to the minimum and warn
        if numel(stimFrameIdx) ~= numel(stimPosVec)
            warning('Recording %d: stimTimes (%d) and stimPosVec (%d) differ in length — trimming to min length.', ...
                    n, numel(stimFrameIdx), numel(stimPosVec));
            L = min(numel(stimFrameIdx), numel(stimPosVec));
            stimFrameIdx = stimFrameIdx(1:L);
            stimPosVec   = stimPosVec(1:L);
        end
    
        % build logical mask of valid trials 
        validTrials = ~isnan(stimFrameIdx) & ~isnan(stimPosVec);
    
        % filter to only valid trials
        stimFrameIdx = stimFrameIdx(validTrials);
        stimPosVec   = stimPosVec(validTrials);
        
        if isempty(stimFrameIdx)
            continue;
        end
        
        % preallocate local trial matrix (trials x window)
        nTrialsLocal = numel(stimFrameIdx);
        localMat = nan(nTrialsLocal, nWindow);
        localOrient = nan(nTrialsLocal,1);
        
        winOffsets = -pre:post;
        for tt = 1:nTrialsLocal
            centerIdx = stimFrameIdx(tt);
            windowIdx = centerIdx + winOffsets;
            
            % clip and fill with NaN where outside bounds
            valid = windowIdx >= 1 & windowIdx <= nFrames;
            tmp = nan(1, nWindow);
            tmp(valid) = trace(windowIdx(valid));
            localMat(tt,:) = tmp;
            
            % map orientation label to index
            pos = stimPosVec(tt);
            if pos == -90
                localOrient(tt) = 1;    % left
            elseif pos == 0
                localOrient(tt) = 2;    % center
            elseif pos == 90
                localOrient(tt) = 3;    % right
            else
                % unknown orientation: add as NaN and store as a separate code
                localOrient(tt) = NaN;
            end
        end
        
        % append to global arrays
        psthData       = [psthData; localMat]; %#ok<AGROW>
        orientationIdx = [orientationIdx; localOrient]; %#ok<AGROW>
        dayIdx         = [dayIdx; repmat(n, size(localMat,1), 1)]; %#ok<AGROW>
        
        % readable labels for trialInfo
        for tt = 1:numel(localOrient)
            if localOrient(tt) == 1, ol = 'left'; 
            elseif localOrient(tt) == 2, ol = 'center';
            elseif localOrient(tt) == 3, ol = 'right';
            else ol = 'other'; end
            orientationLab{end+1,1} = ol; %#ok<SAGROW>
            % store the rec day string if passiveM has it
            if isfield(passiveM, 'day') && numel(passiveM) >= n
                recDayLabels{end+1,1} = passiveM(n).day; %#ok<SAGROW>
            else
                recDayLabels{end+1,1} = sprintf('rec%d', n); %#ok<SAGROW>
            end
        end
        
        % Grab the learning day from aniDates to produce LearnDayIdx
        AL = find(strcmp(animalID, aniDates(:,1)));
        if ~isnan(aniDates{AL,2})
            learnDay = datetime(aniDates{AL,2});
            idxLearn = find({passiveM.day} == learnDay);
            learnDayIdx = dayIdx - idxLearn;
        else 
            learnDayIdx = nan([size(dayIdx)]);
        end

        % Now to make the megatrial matrix
        nTrialsLocal = size(localMat,1);
        
        % append to pooled arrays
        psthData_all       = [psthData_all; localMat];                %#ok<AGROW>
        orientationIdx_all = [orientationIdx_all; localOrient];       %#ok<AGROW>
        mouseIdx_all       = [mouseIdx_all; repmat(m, nTrialsLocal,1)]; %#ok<AGROW>
        
        % learn-day relative index for this recording (n is recording index)
        if exist('idxLearn','var') && ~isempty(idxLearn) && ~isnan(idxLearn)
            rel = n - idxLearn;   
        else
            rel = NaN;            
        end
        learnDayIdx_all = [learnDayIdx_all; repmat(rel, nTrialsLocal,1)]; %#ok<AGROW>
    end
end

[dayOriMeans, groups] = ap.groupfun(@(x) mean(x,1,'omitnan'), psthData_all, [learnDayIdx_all orientationIdx_all]);
[dayOriSem, ~] = ap.groupfun(@(x) nanstd(x, 0, 1) ./ sqrt(sum(~isnan(x),1)), psthData_all, [learnDayIdx_all orientationIdx_all]);

% plotting
relDays = unique(groups(:,1));
nRel = numel(relDays);

figure;
tiledlayout(1, 7, 'TileSpacing','compact','Padding','compact'); % arrange tiles

for di = 8:14
    rel = relDays(di);      
    nexttile;
    hold on;

    % logical for this relative day across groups
    dayMask = groups(:,1) == rel;

    % for each orientation 
    colors = {'b','k','r'};
    labels = {'left','center','right'};
    plotted = false(3,1);

    for ori = 1:3
        mask = dayMask & (groups(:,2) == ori);
        if any(mask)
            M = dayOriMeans(mask, :);  
            SEM = dayOriSem(mask, :);
            % plot filled SEM and mean
            xl = [x, fliplr(x)];
            yl = [M + SEM, fliplr(M - SEM)];
            fill(xl, yl, colors{ori}, 'EdgeColor','none', 'FaceAlpha', 0.2);
            plot(x, M, colors{ori}, 'LineWidth', 1.2);
            plotted(ori) = true;
        end
    end

    xline(0,'--');
    yline(0,'Color','k','Alpha',0.5);
    xlabel('Frames (relative to stimulus)');
    ylabel('Pupil diameter derivative (arb)');
    ylim([-0.05, 0.05]);
    title(sprintf('Learn day = %d', rel));

    hold off;
end

%% Package into a beautiful table

animalID = {};
for i = 1:numel(mouseIdx_all)
    animalID{i,1} = mouseIDs{mouseIdx_all(i)};
end

orientationDir = [];
for y = 1:numel(orientationIdx_all)
    if orientationIdx_all(y,1) == 1
            orientationDir(y,1) = -90;    % left
    elseif orientationIdx_all(y,1) == 2
            orientationDir(y,1) = 0;    % center
    elseif orientationIdx_all(y,1) == 3
            orientationDir(y,1) = 90;     % right
    end
end

pupilDiamByFrame = psthData_all;

pupils = table(animalID, mouseIdx_all, learnDayIdx_all, orientationDir, pupilDiamByFrame);