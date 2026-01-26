%% Data demo
%
% A demo on how to use the server and load data
% Whenever you see a line of code, run it (highlight and press F9) and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% Clone these repos and add to the Matlab path 
% (Home > Set Path > Add with Subfolders...)
%
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab
% https://github.com/kwikteam/npy-matlab

%% File structure and finding recordings

% To connect to the server: 
% - Open up the file explorer, select This PC in the left hand menu 
% - Right click on this PC and select map network drive 
% - Enter the following address: \\qnap-ap001.dpag.ox.ac.uk\APlab\ 
% - It should prompt you for credentials, please enter your MSD details in the following format: 
% - Username – MSD\(MSD username) 
% - Password – normal MSD password 

% The server contains these top folders:
% \Data: for all raw and preprocessed data (all users can write)
% \Users: all users have a folder for their own use (given users can write)
% \Lab: files for the general lab and rigs (only Andy can write)

% Data is organized on the server in these folders: 
% \Data
% |- animal name (initials & number, e.g. AP012)
% | |- recording date (yyyy-mm-dd, e.g. 2023-03-16)
% |  |- recording time (Recording_HHMM, e.g. Recording_1021 started at 10:21am)
% |  | |- recording data by modality (e.g. mousecam, widefield)
% |  |- day data spanning multiple recordings (e.g. ephys, widefield)
% |- non-recording data (e.g. histology)


%% Finding data from recordings

% The function 'plab.find_recordings' is used to find data from recordings
% The syntax is: plab.find_recordings(animal,recording_day,workflow)
% (type 'help plab.find_recordings' for documentation).
% 'animal' must be filled in, other specifics can be left out to return all
% relevant recordings.

% --- Working with plab.find_recordings

% This line finds all recordings done with animal AP010:
% (day and workflow are not entered)
recordings = plab.find_recordings('AP010');

% 'recordings' is a structure with one entry for each day containing:
% - day: day of recording
% - index: index of recording within the day for that animal (e.g. 3rd recording = [3])
% - recording: time of recording (HHMM)
% - workflow: Bonsai workflow (i.e. stimuli or task protocol)
% - mousecam: if there was mouse camera recording (0/1 = no/yes)
% - widefield: if there was widefield recording (0/1 = no/yes)
% - ephys: if there was electrophysiology recording (0/1 = no/yes)

% [EXERCISE]
% 1) How many days were recordings performed on AP010?
days = {recordings.day}; % Grab array of days
RecordingDays = numel(days); % Assign RecordingDays to number of elements
disp(RecordingDays); % Tell us
% Ans: 16

% 2) How many recordings were performed on AP010 across all days?
indexes = {recordings.index}; % Grab array of indexes

totalRecordings = 0; % Initialize empty array of total recordings

for i = 1:numel(indexes) % For each index entry...
    totalRecordings = totalRecordings + numel(indexes{i}); % Add number of elements to count
end
disp(totalRecordings) % Tell us

% Notes from 16/10/2025: this loop, where totalRecordings is assigned
% values in separate places, is bad practice bc it will cause mismatches when
% troubleshooting. That's why the below is better. 

% OR

totalRecordingsAlt = sum(arrayfun(@(x) numel(x.index), recordings)); % Uses an array function to ask for the number of elements in each index entry from recordings, then sums them
disp(totalRecordingsAlt); % Tell us



% Ans: 37

% 3) How many unique Bonsai workflows were used, and what were they?

allWorkflows = {}; % Initialize list

for i = 1:numel(recordings)
    allWorkflows = [allWorkflows; recordings(i).workflow(:)]; % Add each workflow in a clunky list
end

uniqueWorkflows = unique(allWorkflows);
disp("There were "+numel(uniqueWorkflows)+" workflows:");
disp(uniqueWorkflows)

% OR

allWorkflowsAlt = vertcat(recordings.workflow) % Vertically concatenate all of the workflow entries

uniqueWorkflowsAlt = unique(allWorkflowsAlt) % Get unique workflows

disp("There were "+numel(uniqueWorkflowsAlt)+" workflows:"); % Tell us
disp(uniqueWorkflowsAlt)

fprintf("There were %d workflows", numel(uniqueWorkflowsAlt))
fprintf("Workflows = %s\n", uniqueWorkflowsAlt{:})

%fprintf is more efficient and cool. The % operators differ in what they
%accept, and \n starts a new line

% Ans: 5, listed above

% This line finds all recordings of AP010 on 2023-08-24: 
% (workflow is not entered)
recordings = plab.find_recordings('AP010','2023-08-23');

% [EXERCISE]
% Which recording modalities were recorded this day?

modalities = {'mousecam','widefield','ephys'}; % List what fields we're asking about

areModalitiesRecorded = cellfun(@(f) any(recordings.(f) ~= 0), modalities); % We only need to know if a given field isn't a 0

% This format any(recordings.(f) ~= 0) is redundant; the section inside the
% bracket would be done by 'any'

disp(modalities(areModalitiesRecorded)); % Tell us which came back as recorded on this day

% Ans: All three. Tested on differing days and didn't break. 

% This line finds all recordings of AP010 with workflow 'lcr_passive':
% (date is not entered)
recordings = plab.find_recordings('AP010',[],'lcr_passive');

% [EXERCISE]
% 1) How many dayss of this workflow included widefield? 

% Since we've already filtered, so to speak, straightforward? 
totalWidefieldRecordings = sum(cat(2, recordings.widefield))
% Ans: 17

% Alternate for how many days had widefield
daysWidefield = sum(cellfun(@any, {recordings.widefield}))

% 2) How many recordings of this workflow included electrophysiology? 

totalEphysRecordings = sum(cellfun(@(e) sum(e(:)), {recordings.ephys})); % same again
% Ans: 3

% 3) Make a new variable which is a subset of 'recordings' with widefield

hasWidefield = cellfun(@(w) any(w(:)), {recordings.widefield}); % Seeing which has any widefield
widefieldRecordings = recordings(hasWidefield); % Subset

% Not exactly sure if this fits the bill - its more like a subset of days
% which contain at least one widefield recording. 

% This line finds recordings of AP010 with workflow 'lcr_passive' on
% 2023-08-16 (all fields are entered)
recordings = plab.find_recordings('AP010','2023-08-16','lcr_passive');

% [EXERCISE]
% This day has 2 recordings of the same workflow. This usually is
% because there was an issue with the first one and it needed to be re-run.
% In the number order of that day's recordings, which ones were
% 'lcr_passive'? (hint: 'index')
disp(recordings.index);
% Ans: the first and the third


% The 'workflow' can include * as a wildcard. For example, in the task
% workflow 'stim_wheel_right', there is a 'stage1' and 'stage2' version.
% This returns specifically recordings with stage1:
recordings = plab.find_recordings('AP010',[],'stim_wheel_right_stage1');
% And this returns recordings with either stage: 
recordings = plab.find_recordings('AP010',[],'stim_wheel_right*');

% [EXERCISE]
% Write code to return the day AP010 switched from stage1 to stage2.

% This phrasing has me scared I wasn't supposed to be writing code
% before...

% Assuming they're always sequential, easiest way is just to grab the first day
% where stage2 appears. 
recordings = plab.find_recordings('AP010', [], 'stim_wheel_right_stage2');

disp(recordings(1).day);

% Exercise to do it the hard way

firstDay = cellfun(@(g) g(end) == '2', cat(1, recordings.workflow), 'UniformOutput', true)

firstDay = cellfun(@(g) all(g == 'stim_wheel_right_stage2'), cat(1, recordings.workflow))

recordings(find(firstDay, 1)).day

% Ans: 2023-08-13

%% Constructing server paths

% The class 'plab.locations' is used to find or construct paths on the server.

% --- Standardized generic locations

% This is the general path to the server:
plab.locations.server_path
% This is the data path on the server:
plab.locations.server_data_path
% This is where local data should be kept on each computer: 
plab.locations.local_data_path

% [EXERCISE]
% Use the function 'fullfile' and the above information to construct the
% path: server/Users/(your name)/test_path

serverPath = plab.locations.server_path; % Grab the server path

myName = "Peter_Gorman"; % My name as on server

fullPath = fullfile(serverPath, 'Users', myName, 'test_path'); % Construct

disp(fullPath)

% --- Recording locations

% The method 'plab.locations.filename' is used to construct paths/filenames
% for recording data
% (type 'help plab.locations.filename' for documentation). 
% The syntax is:   constructed_filename = filename('server' or 'local',animal,rec_day,rec_time,folder1,...,folderN)

% For example: 
% This constructs the path for animal AP001 on day 2000-01-01 at 12:00:
plab.locations.filename('server','AP001','2000-01-01','1200')
% Input arguments can be added or removed as necessary, for example:
% This constructs the path for the above recording and the subfolder
% 'mousecam'
plab.locations.filename('server','AP001','2000-01-01','1200','mousecam')
% This constructs the path for the above animal and day in the folder
% 'ephys' (not in a recording time folder):
plab.locations.filename('server','AP001','2000-01-01',[],'ephys')
% And this example constructs a path 'histology' in the animal folder (not
% in a day folder): 
plab.locations.filename('server','AP001',[],[],'histology')

% Note that paths with plab.locations.filename are constructed whether or
% not the folder/file exists (e.g. above example paths do not exist).

% [EXERCISE] 
% Use 'plab.find_recording' to find the recording time for AP010 on
% 2023-08-10 with workflow 'lcr_passive', then use
% 'plab.locations.filename' to construct the path to the 'widefield'
% subfolder within the folder for that recording. Write a line to check
% whether that path exists on the server.

recording = plab.find_recordings('AP010', '2023-08-10', 'lcr_passive');
recordingTime = recording.recording{1}; % Extract the recording time
widefieldPath = plab.locations.filename('server', 'AP010', '2023-08-10', recordingTime, 'widefield');
disp(isfolder(widefieldPath)); 

%% Loading data, and data types

% There isn't standardized lab code for loading data, since this is often
% customized for each person depending on their needs. This block demos my
% code, which can be used as a template, or as-is if it works for you (note
% that it's subject to regular changes).

% My loading code is in my repository (petersaj/AP_scripts_peterslab) and
% is called with 'ap.load_recording' after defining animal/day/time. This
% loads an example dataset with passive stimuli, widefield, and ephys. Run
% this data, then each data type will be explained below:
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
verbose = true; % this turns on/off progress display in command line
ap.load_recording;

% --- Timelite

% "Timelite" is our GUI for recording analog signals with a DAQ. These
% include things like wheel movement, stimulus screen photodiode, camera
% frame captures, etc. It is saved as:
% server\Data\animal\day\time\timelite.mat
% as the structure 'timelite':
timelite

% The structure 'timelite' contains these structures: 
% - daq_info: information about the DAQ recording
% - data: the recorded data as N timepoints x M signals
% - timestamps: the timestamps for each data point as N timepoints x 1

% The names of the recorded signals are under
% 'timelite.daq_info(channel number).channel_name'. For example, the
% photodiode is recorded under channel 5, and the name can be found by:
photodiode_channel = 5;
timelite.daq_info(photodiode_channel).channel_name

% The corresponding data is in the 5th column of 'data'. For example, this
% plots photodiode data by time: 
figure;plot(timelite.data(:,photodiode_channel));
xlabel('Time (s)');
ylabel('Voltage');
title('Photodiode signal');

% (The photodiode is over a square on the screen that turns white/black to
% indicate changes to the stimulus).

% [EXERCISE] 
% 1) Write code to identify the channel number for 'widefield_camera'
% (widefield camera exposures) and plot the data for that channel with
% timestamps on the x-axis

names = {timelite.daq_info.channel_name}; % Getting the channels subfield out first
widefield_channel = find(strcmp(names, 'widefield_camera')); % Correctly tells us it's three

% Plot it by time
figure;plot(timelite.timestamps, timelite.data(:,widefield_channel));
xlabel('Time (s)');
ylabel('Signal');
title('Widefield signal');

% This seems... not quite useful

% 2) Find the timestamps when the widefield camera starts each exposure
% (first sample when the exposure signal is high after being low)

% This signal quality is very consistent so I'm tempted to just manually
% threshold 

widefieldSignal = timelite.data(:, widefield_channel);

isHigh = widefieldSignal > 1.5; % Seeing when signal is arbitrarily high

risingIdx = find(diff([double(isHigh)]) == 1) + 1 ;   % Indices of first high sample per burst
risingTimes = timelite.timestamps(risingIdx);

% 3) Check that the timestamps in 2 match the variable
% 'widefield_expose_times' (this is created in ap.load_recording)

diffs = risingTimes - widefield_expose_times; % Doing this by checking the difference between the two arrays. 
maxDiff = max(abs(diffs));
fprintf('Max difference: %.12g s\n', maxDiff);

isequal(risingTimes, widefield_expose_times)

% Close to a match...

% This made me curious about how the arbitrary threshold impacts the
% accuracy of my approach. Here's code to try different thresholds and see
% which is best. 

timestamps = timelite.timestamps(:);

thresholds = 0.05:0.05:2.95;   % Range of thresholds to test
maxDiffs = nan(size(thresholds));  % This approach could be improved

for i = 1:numel(thresholds)
    th = thresholds(i);

    % Detect rising edges at this threshold
    isHigh = widefieldSignal > th;
    risingIdx = find(diff([0; double(isHigh)]) == 1);
    risingTimes = timestamps(risingIdx);

    % Compare lengths as long as both have values
    if isempty(risingTimes) || isempty(widefield_expose_times)
        maxDiffs(i) = NaN;
        continue
    end

    % Handle unequal lengths
    n = min(numel(risingTimes), numel(widefield_expose_times));
    diffs = risingTimes(1:n) - widefield_expose_times(1:n);
    maxDiffs(i) = max(abs(diffs));
end

% Find the smallest max difference
[bestDiff, bestIdx] = min(maxDiffs);
bestThreshold = thresholds(bestIdx);

fprintf('Best threshold: %.2f (maxDiff = %.6f s)\n', bestThreshold, bestDiff);

% Plot results
figure;plot(thresholds, maxDiffs, 'o-');
xlabel('Threshold value');
ylabel('Max difference (s)');
title('Thresholds for best rising edge alignment');
grid on;

% Looks like 1.9 - 2.3 is the way to go. Cool

% One important signal recorded in Timelite is the wheel which the mouse
% can turn left and right. This is in the channel called 'wheel', which
% represents the position of the wheel relative to the start of the
% recording. Besides position, it is also helpful to have wheel velocity,
% and binary classification of whether the wheel is moving or not. This
% information is calculated by 'ap.parse_wheel', and outputs these
% variables - 
% wheel_velocity: velocity of the wheel for each Timelite timestamp
% wheel_move: binary vector representing movement (1) or quiescence (0)

% [EXERCISE] 
% 1) Plot time vs. raw Timelite wheel position and the wheel velocity on
% separate axes. Using 'wheel_move', plot the wheel velocity only when the
% wheel is moving in a separate color (by setting quiescent times to NaN,
% which are not plotted).

wheel_channel = find(strcmp(names, 'wheel')); % Have to grab channel number first

% Create mask for changing velocity colour
wheel_velocity_move = wheel_velocity;
wheel_velocity_move(~wheel_move) = NaN; % NaN out 

% Plot all together
figure;
subplot(2,1,1);
plot(timestamps, timelite.data(:, wheel_channel), Color = 'k');
title('Raw wheel position');
subplot(2,1,2);
hold on;
plot(timestamps, wheel_velocity, Color=[.7 .7 .7]);
plot(timestamps, wheel_velocity_move, Color="r");
hold off;
title('Wheel velocity');

% 2) Using 'wheel_move', find the onset and offset times of all movements.
% Plot these as lines (colored differently for onset and offset) on your
% previous velocity plot.

d = diff(wheel_move);
onIdx  = find(d == 1) + 1;         % Indices of first sample of each movement 
offIdx = find(d == -1);    % Indices of last sample of each movement

% Convert to times
onTimes  = timestamps(onIdx);
offTimes = timestamps(offIdx);

% Just the velocity plot
figure;
hold on;
plot(timestamps, wheel_velocity, Color=[.7 .7 .7]);
plot(timestamps, wheel_velocity_move, Color="r");
xline(onTimes, Color='g', LineWidth=0.1, LineStyle='--');
xline(offTimes, Color='m', LineWidth=0.1, LineStyle='--');
hold off;
title('Wheel velocity');

% Truly offensive to the eyes

% --- Bonsai

% We use the program Bonsai to run our stimuli and tasks
% (https://bonsai-rx.org/). Bonsai is an actively supported, highly
% customizable, visual interfaced framework for designing experiments. 

% Bonsai files are called 'workflows'. ap.load_recording loads the name of
% the currently loaded workflow as 'bonsai_workflow'. This one is
% 'lcr_passive', which is passive stimuli presented on the left, center,
% and right screens.
bonsai_workflow

% Bonsai can change how it saves data, and loading scripts can specify how
% that data is loaded and parsed. In this demo, ap.load_recording loads
% data from Bonsai as 'trial_events':
trial_events
% which is a structure containing:
% - parameters: parameter values corresponding to the whole workflow, e.g.
% the duration that each stimulus was displayed was:
trial_events.parameters.StimDuration
% - values: an N trials x 1 array of values of saved events, e.g. the
% order of X-positions for the 3 stimuli presented on trial 4 was:
trial_events.values(4).TrialStimX
% - timestamps: an N trials x 1 array of timestamps of saved events from
% 'values', e.g. the time each of 3 stimuli was turned on/off on trial 4
% was:
trial_events.timestamps(4).StimOn

% Note that timing information in Bonsai is only approximate (when the
% software gave the command), not actual (when the stimulus was physically
% drawn on the screen). All time information should be used from Timelite,
% and only qualitative information should be used from Bonsai (e.g. which
% stimulus was presented on a trial).

% [EXERCISE] 
% The stimulus onset times are loaded by ap.load_recording as
% 'stimOn_times'. Write code to pull out a subset of these timestamps
% corresponding to a stimulus X position (TrialStimX) of 90 (meaning it was
% on the right-hand screen).

stimPositionsVec = vertcat(trial_events.values.TrialStimX); % Grab the cells of all stimulus positions

idx90 = find(stimPositionsVec == 90); % Find identities for each position in that vector which is 90

rightHandPresentations = stimOn_times(idx90); % Apply this order to the timelite array

% disp(rightHandPresentations);

% --- Mousecam

% We record video of the front of the mouse during all experiments. The
% filename is loaded in ap.load_recording as:
mousecam_fn
% and the timestamps of each frame (in Timelite clock) is:
mousecam_times

% Mousecam images can be read into Matlab with a VideoReader object. For
% example:
mousecam_vr = VideoReader(mousecam_fn); % create VideoReader object

load_frame = 1; % define frame to read
mousecam_im = read(mousecam_vr,load_frame); % read frame

figure; % create a figure
imagesc(mousecam_im); % draw the image (sc = scaled colors)
axis image % make the x/y axes have equal aspect ratios
colormap('gray'); % set the colormap to gray

% You can also load in multiple frames at the same time by defining the
% start/end frame. For example:
load_frames = [1,10]; % define frame interval to read
mousecam_im = read(mousecam_vr,load_frames); % read frame
% (the VideoReader by default loads in multiple images in the 4th
% dimension, allowing for colors in the 3rd dimension. Our images are
% grayscale, so we can use 'squeeze' to remove the singleton 3rd dimension.
% Look at the size of the natively loaded data: 
size(mousecam_im)
% and then the size of the 'squeezed' data:)
mousecam_im = squeeze(mousecam_im); 
size(mousecam_im)
% This is an example of how to view a 3D matrix using my function
% ap.imscroll, which plots each 'page' (dim 1+2) which can be scrolled
% through in dim 3:
ap.imscroll(mousecam_im);

% [EXERCISE] 
% Create an average mousecam movie for -20:+20 frames around stimulus X =
% 90 presentations, and separately for stimulus X = -90. Is there a
% difference in behavior when these two stimuli are presented?

% First lets get the presentations. We have right-hand already.

RightPres = rightHandPresentations;
LeftPres = stimOn_times(find(stimPositionsVec == -90)); 

% We can use these with mousecam_times but we need to define a -20:+20
% window for each by approximating time to frame number. They don't quite
% match
rightFrameIdx = interp1(mousecam_times, 1:numel(mousecam_times), RightPres, 'previous', NaN);
leftFrameIdx = interp1(mousecam_times, 1:numel(mousecam_times), LeftPres, 'previous', NaN);

% Now we create ranges for each presentation

rightRange = rightFrameIdx + (-20:20); 
leftRange = leftFrameIdx + (-20:20); 

% Averaging the frame images across each *column* in the matrix rightRange

vid = VideoReader(mousecam_fn); % Giving it my own name to not get mixed up

[nPres, nWindow] = size(rightRange); % Need the numbers of presentations and window size

firstframe = [size(read(vid,1)), nWindow];

avgRight = zeros(firstframe, "double"); % Initialising. Not ideal :/

for i = 1:nWindow
    frameIds = [rightRange(:,i)]; % Get the frame indices for column in question
    
    sumFrames = zeros(448, 496, 'double'); % Cheating slightly bc we already know the dimensions from above

    for j = 1:numel(frameIds) % Stack em
        frameIdx = frameIds(j);
        frame = read(vid, frameIdx);
        sumFrames = sumFrames + double(frame);
    end

    avgRight(:,:,i) = sumFrames /numel(frameIds); % Average pixel intensities into the avgRight matrix
    % This strikes me as a strange way to do this but if it works it works
end

averageRight = zeros(firstframe, "double");

for i = 1:nPres
    divPres = double(squeeze(read(vid, rightRange(i, [1, end]))))/nPres;
    averageRight = averageRight + double(divPres);
    disp(i);
end



% Okay well now it's not a movie anymore. Go back to uint8 format?

avgRight_uint8 = uint8(averageRight);
implay(avgRight_uint8); % Didn't expect either of these to actually work

ap.imscroll(averageRight);

% % How I had it before
% currAxes = axes;
% 
% for k = 1:size(avgRight,3)
%     image(avgRight(:,:,k),"Parent",currAxes)
%     currAxes.Visible = "off";
%     colormap gray
%     pause(1/v.FrameRate)
% end

% Okay, now I'll do it all again for the left presentations. Horribly
% clunky but maybe we can refactor

avgLeft = zeros(448, 496, nWindow, "double"); % Initialising left. Not ideal :/

for i = 1:nWindow % Technically this is the nWindow for the rightRange again.
    frameIds = [leftRange(:,i)]; % Get the frame indices for column in question
    
    sumFrames = zeros(448, 496, 'double'); % Cheating slightly but we already know the dimensions from above

    for j = 1:numel(frameIds) % Stack em
        frameIdx = frameIds(j);
        frame = read(vid, frameIdx);
        sumFrames = sumFrames + double(frame);
    end

    avgLeft(:,:,i) = sumFrames /numel(frameIds); % Average pixel intensities into the avgRight matrix
end

avgLeft_uint8 = uint8(avgLeft);
implay(avgLeft_uint8);

leftminusright = avgLeft_uint8 - avgRight_uint8;

max(leftminusright, [], "all") % So we can adjust max values to see difference better

implay(leftminusright); % This doesn't tell me a thing. Anyways

% To answer the question, "is there a difference in behavior when these two
% stimuli are presented", we'd need to get some idea of
% presentation-presentation variability (better than just an average) and
% get some model which identifies features better

% --- Widefield and electrophysiology

% If available in a recording, ap.load_recording will load and prepare
% widefield and electrophysiology data. For demos on working with that
% data, see: 
% - plab_widefield_demo
% - plab_ephys_demo

% I have a function to scroll through experiment data for exploratory
% purposes, which displays the mousecam (left), widefield (center), and
% multiunit ephys (right), which can be scrolled through time:
ap.expscroll

% [EXERCISE] 
% There is a drop of sucrose available to the mouse at the beginning of
% this recording (left over from the task done previously). Which part of
% the ephys probe has increased activity when the mouse consumes this?

% By inspection, I'm saying the 2500-3200 range

% --- Loading partial datasets

% Recordings can be part-loaded by ap.load_recording if only some data is
% necessary by creating a 'load_parts' structure, which can toggle loading
% with fields 'widefield','ephys','mousecam'. Timelite and behavior are
% always loaded. If 'load_parts' exists, any field not defined will default
% to not loading. For example, in this dataset, you can turn off loading of
% widefield and ephys data by doing:
clear all % (clears data loaded above)
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
load_parts.mousecam = true; % (this is the new line)
verbose = true;
ap.load_recording; % (this is much faster now and omits widefield and ephys)




















