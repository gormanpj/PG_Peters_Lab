
animal = 'PG003';
rec_day = '2026-03-25';
rec_time = '1741';
verbose = true;
ap.load_recording;

stim_x = vertcat(trial_events.values.TrialStimX);
ap.cellraster(stimOn_times,stim_x);


stim_y = vertcat(trial_events.values.TrialStimY);

ap.cellraster(stimOn_times,stim_y);

stim_type = vertcat(trial_events.values.TaskType);
stim_type = stim_type(1:numel(stimOn_times));

ap.cellraster(stimOn_times,stim_type);
