% Exploratory behaviour script

% Visualize averaged wheel movement around stimulus
animal = 'PG001'; 
use_workflows = {'stim_wheel_right_stage\d_noincorrect','stim_wheel_right_stage2'};
recordings = [];
for i = 1:numel(use_workflows)
    rec = plab.find_recordings(animal,[],use_workflows{i});
    rec = reshape(rec, 1, []);
    recordings = [recordings, rec]; %ok<AGROW>
end

x = (-5000:5000);

figure; hold on;
title(animal, use_workflows, 'Interpreter', 'none');
axis square;
colours = {[.3 0 .7], [.5 0 .7], [.5 .2 .7], [.5 .4 .7], [.5 .6 .7], [.5 .6 .9],[.5 .7 1], [.5 .9 1], [.6 1 1],[.8 .8 1], [1 1 1]};
for v = 1:length(recordings);
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
    
    wheelStim = mean(wheel_move(stimRanges), 1);

    % Percentage of time with stim on
    stimTotal = 100*(sum(photodiode_trace >= 3) / length(photodiode_trace));

    plot(x, wheelStim, 'Color', colours{v} ,'LineWidth', 3);
    xticks(-5000:1000:5000)
    xticklabels(-5:1:5)
    ylabel('Wheel movement %');
    xlabel('time around stimulus (s)');
    xline(0, '--');
    ylim([0 1]);
    text(5500,0.95 - 0.05*v, [rec_day ' (' num2str(round(stimTotal, 3)) ' Stim On %)'], "FontSize",14, "Color", colours{v});
end
hold off;

% Plot wheel movements relative to stimuli 


figure; tiledlayout(length(recordings),2,'TileIndexing', 'columnmajor');
title(animal, use_workflow);
for v = 1:length(recordings);
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
    
    wheelPos = cumtrapz(wheel_velocity);
    nexttile; 
    plot(wheel_velocity, 'Color', 'k');
    xline(stimFrameIdx, 'Color','r');

    nexttile;
    plot(wheelPos, 'k');
end
