% Exploratory behaviour script

% Visualize averaged wheel movement around stimulus
animal = 'PG001'; 
use_workflows = {'stim_wheel_right_AVAP_noincorrect'};
recordings = [];
for i = 1:numel(use_workflows)
    rec = plab.find_recordings(animal,[],use_workflows{i});
    rec = reshape(rec, 1, []);
    recordings = [recordings, rec]; %ok<AGROW>
end

x = (-5000:5000);
colours = {[.3 0 .7], [.5 0 .7], [.5 .2 .7], [.5 .4 .7], [.5 .6 .7], [.5 .6 .9],[.5 .7 1], [.5 .9 1], [.6 1 1],[.8 .8 1], [1 1 .8], [1 1 .6], [1 1 .4], [1 1 .2]};

figure;

% 1 row, 3 columns
axAgg = subplot(1,3,1); hold(axAgg,'on'); 
title(axAgg,'Aggregate'); 
ylabel(axAgg,'Wheel movement %'); 
xline(axAgg,0,'--'); 
ylim(axAgg,[0 1]); 
xticks(axAgg,-5000:1000:5000); 
xticklabels(axAgg,-5:1:5);
axis(axAgg,'square');

axAP  = subplot(1,3,2); hold(axAP,'on'); 
title(axAP,'Appetitive'); 
xline(axAP,0,'--'); 
ylim(axAP,[0 1]); 
xticks(axAP,-5000:1000:5000); 
xticklabels(axAP,-5:1:5);
axis(axAP,'square');

axAV  = subplot(1,3,3); hold(axAV,'on'); 
title(axAV,'Aversive'); 
xline(axAV,0,'--'); 
ylim(axAV,[0 1]); 
xlim(axAV,[-5000 6500]); 
xticks(axAV,-5000:1000:5000); 
xticklabels(axAV,-5:1:5);
axis(axAV,'square');

xlabel(axAgg,'time (s)');
xlabel(axAP,'time (s)');
xlabel(axAV,'time (s)');

for v = 2:length(recordings)
    rec_day = recordings(v).day;
    rec_time = recordings(v).recording{end};
    verbose = false;
    load_parts = struct;
    load_parts.behavior = true;
    ap.load_recording

    odds = mod(1:length(photodiode_times), 2);
    evens = mod(2:length(photodiode_times)+1, 2);

    stimOnTimes = photodiode_times(odds==1);
    stimOffTimes = photodiode_times(evens==1);

    stimFrameIdx = interp1(timelite.timestamps, 1:numel(timelite.timestamps), stimOnTimes, 'previous', NaN);

    stimRanges = stimFrameIdx + (-5000:5000);
    [nStim, nWindow] = size(stimRanges);

    % Build per-trial wheel trace matrix (no guards)
    wheelPerTrial = nan(nStim, nWindow);
    for t = 1:nStim
        wheelPerTrial(t, :) = wheel_move(stimRanges(t, :));    % assume indices valid
    end

    % Split to aversive/appetitive
    TaskType = vertcat(trial_events.values.TaskType);
    TaskType = TaskType(1:numel(stimOffTimes));

    % Means 
    wheelStim = mean(wheelPerTrial, 1);
    wheelStimAP = mean(wheelPerTrial(TaskType==0, :), 1);
    wheelStimAV = mean(wheelPerTrial(TaskType==1, :), 1);

    % Percentage of time with stim on
    stimTotal = 100*(sum(photodiode_trace >= 3) / length(photodiode_trace));

    % Plot into the three subplots (same colour per recording)
    plot(axAgg, x, wheelStim, 'Color', colours{v}, 'LineWidth', 2);
    plot(axAP,  x, wheelStimAP, 'Color', colours{v}, 'LineWidth', 2);
    plot(axAV,  x, wheelStimAV, 'Color', colours{v}, 'LineWidth', 2);

    % Labels
    text(axAV, 5200, 0.95 - 0.05*v, ...
    [rec_day ' (' num2str(round(stimTotal, 3)) ' % On)'], ...
    "FontSize",12, "Color", colours{v}, ...
    "HorizontalAlignment","left");
end
% 
% Tgoodbad = vertcat(trial_events.timestamps.TaskType)
% bad = Tgoodbad(vertcat(trial_events.values.TaskType) == 1);
% good = Tgoodbad(vertcat(trial_events.values.TaskType) == 0);
% 
% figure;hold on
% plot(photodiode_bw_interp, 'k')
% 
% plot(timelite.data(:,7), 'b')
% xline((seconds(good-start)*1000)+14400, 'g')
% xline((seconds(bad-start)*1000)+14400, 'r')