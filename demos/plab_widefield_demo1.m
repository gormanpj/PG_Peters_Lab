%% Widefield demo 1
%
% A demo on how to work with widefield data
% (do plab_data_demo first)
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% These repositories are necessary for this demo:
%
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab
% https://github.com/kwikteam/npy-matlab

%% Example dataset

% My loading script (ap.load_recording) loads and preprocesses widefield
% data (specifically in ap.load_widefield), so that will be used for this
% demo.

% Load this data for the following 2 sections: 
animal = 'AP010';
rec_day = '2023-08-07';
rec_time = '1616';
verbose = true;
ap.load_recording;

%% Working with widefield data

% Widefield data is captured in pixels at ~400x400, which is 160k pixels
% per frame. This can be hard to work with, because it takes up a lot of
% memory and operations across this many pixels can be slow or impossible.
% For example, the loaded data is just 5 minutes of imaging, and the raw
% widefield images are 8.6GB. The raw values are also 16-bit numbers, but
% many matlab operations require 32-bit, so that requires 8.6GB * 2 =
% 17.2GB RAM just for this short dataset.

% To make the data more accessible, we compress it with singular value
% decomposition (SVD). Instead of each pixel being independent, SVD
% expresses the data as weighted maps of all pixels (U, spatial components)
% which vary over time (V spatial components). The data can then be
% reconstructed by multiplying corresponding U's and V's and summing across
% components, which can be done with matrix multiplication data = U*V.

% The first component explains the most variance in the data, and
% subsequent components explain progressively less variance. This means we
% can either use all components (= number of pixels or timepoints,
% whichever's smaller) to exactly reconstruct the data, or a subset of
% components to reconstruct most of the variance. We keep 2000 components,
% since this explains >95% of the variance. This brings the memory in our
% example from 8.6GB to 1.5GB.

% --- Working with widefield SVD data

% ap.load_recording loads these components as: 
% wf_U: U (spatial components - Y pixels x X pixels x N components)
% wf_V: V (temporal components - N components x M timepoints)
% wf_t: timestamps in seconds for each frame (M timepoints x 1, in Timelite clock)

% Take a look at the U's (spatial components). Note the early components
% have structure (patterns of pixels that are correlated/anti-correlated
% and explain much of the variance - the first component is the average by
% definition), and later components look like noise (these pick up small
% variance):
ap.imscroll(wf_U);
clim([-0.01,0.01]); 
axis image
colormap(gray);

% Each page of U is a component, which corresponds to the row in V. For
% example, this plots the U and V for component 3:
plot_component = 3;

figure; colormap(gray)
tiledlayout(1,2,'tilespacing','tight');

nexttile;
imagesc(wf_U(:,:,plot_component));
axis image
title(sprintf('U component %d',plot_component))
nexttile; 
plot(wf_t,wf_V(plot_component,:));
title(sprintf('V component %d',plot_component))

% The U's (spatial) and V's (temporal) are discussed above, but SVD
% produces a third matrix S (amplitude of each component), such that
% standard SVD expresses data D as D = USV. For simplicity, we save V's as
% S*V. Because of this, our early V components have larger amplitudes than
% late components because they drive more variance. Here's an example of
% amplitude differences between component 3 and 500, note that 3 is larger
% than 500:
plot_components = [3,500];
figure;
plot(wf_t,wf_V(plot_components,:));
xlabel('Time (s)');
title(sprintf('V components: %d,%d',plot_components(1),plot_components(2)));

% Pixel values for each frame can be reconstructed by matrix multiplication
% of U's and V's (this requires some reshapes from how we normally store
% them: we store U as Y x X pixels x components, but for matrix
% multiplcation it has to be flattened into pixels x components). Here's an
% example of reconstructing one frame:
use_frame = 500;
wf_U_flat = reshape(wf_U,[],size(wf_U,3)); % flatten U for U*V operation
example_fluorescence_flat = wf_U_flat*wf_V(:,use_frame); % matrix multiply all components for given frame
example_fluorescence_frame = reshape(example_fluorescence_flat,size(wf_U,[1,2])); % reshape flat pixels into frame size
figure;
imagesc(example_fluorescence_frame);
axis image off;
title('Example frame fluorescence');

% Note that the image above looks different from the raw image, because it
% is mean-subtracted. The average raw image is also saved and loaded, and
% can be displayed like this:
figure;
imagesc(wf_avg);
axis image off;
title('Average widefield image');

% As mentioned above: the benefit of SVD is that you don't need to use all
% components, because you can capture most of the variance with a subset.
% For example, here we can reconstruct the same frame above, but with
% different numbers of components. Note that more components give the image
% more focused regions of activity, but there isn't much difference in the
% images after you have ~50 components:

% (set the numbers of components to use for each reconstruction)
use_components = [5,10,50,100,200,500];
figure;
tiledlayout(1,length(use_components));
for curr_reconstruct = 1:length(use_components)
    % (set the number of components for this reconstruction)
    curr_components = use_components(curr_reconstruct);
    % (this is as above, but indexing the columns of wf_U_flat - the rows
    % are pixel, and the rows of wf_V - the columns are frame)
    example_fluorescence_flat = wf_U_flat(:,1:curr_components)*wf_V(1:curr_components,use_frame);
    example_fluorescence_frame = reshape(example_fluorescence_flat,size(wf_U,[1,2]));

    nexttile;
    imagesc(example_fluorescence_frame);
    clim([-0.006,0.006])
    axis image off
    title(sprintf('%d components',curr_components));
end

% To avoid writing out the matrix multiplication each time, we have the
% function 'plab.wf.svd2px', which does the reshape and multiplication:
use_frame = 500;
example_fluorescence_frame = plab.wf.svd2px(wf_U,wf_V(:,use_frame));
figure;
imagesc(example_fluorescence_frame);
axis image off;
title('Example frame fluorescence');

% Multiple frames can be reconstructed together (or all of them, if you
% don't index the V's). Scroll through this data to see a movie of the
% reconstructed widefield frames:
use_frames = 300:400;
example_fluorescence_frames = plab.wf.svd2px(wf_U,wf_V(:,use_frames));
ap.imscroll(example_fluorescence_frames);
axis image;

% --- Operations on widefield SVD temporal components (V's)

% Another benefit of SVD is that we can do any linear operations (+,-,/,*)
% on the V temporal components, rather than reconstructed pixels, because
% SVD reconstruction is a linear operation. For example, if we wanted to
% average frames together, we can average the V's first and then
% reconstruct, which is equivalent to reconstructing pixels and then
% averaging: 
frames_to_average = 1:10:100; % select frames to average
avg_V = mean(wf_V(:,frames_to_average),2); % average V selected frames for all components
avg_px = plab.wf.svd2px(wf_U,avg_V);
figure;imagesc(avg_px)
axis image off

% [EXERCISE] 
% 1) Show that the above point is true: average frames by (1) averaging
% V's, then reconstructing to pixels (2) reconstructing to pixels, then averaging.
% Write code to check if 1 and 2 are equivalent.

framesToAverage = 2:20:200; % select new frames to average

% First replicate the above but with new interval (1)

avgV = mean(wf_V(:,framesToAverage),2); % average V selected frames for all components
avgPx1 = plab.wf.svd2px(wf_U,avgV);

% Then do the same but with wf_U pixel averages

reconstructedFrames = plab.wf.svd2px(wf_U,wf_V);
avgPx2 = mean(reconstructedFrames(:,:,framesToAverage),3);

% Then plot and use isequal

figure; tiledlayout(1,2);
nexttile;
imagesc(avgPx1);
title('(1) Avg. V, then reconstruct pixels');
nexttile;
imagesc(avgPx2);
title('(2) Reconstruct pixels, then average');
% They look indistinguishable 

isequal(avgPx1, avgPx2); % this actually returns 0 - interesting. 
d = max(abs(avgPx1 - avgPx2), [], 'all');
fprintf('Max difference is %.12f. \n', d);
% So probably in the way things are calculated?


% 2) Show that the order of operations does matter for non-linear
% operations: do the same as above but with standard devation (which is
% non- linear because it uses a square root)

framesToAverageSD = 3:30:300; % new frames to average again

% (1)
sdV = std(wf_V(:,framesToAverageSD), 0 ,2); 
sdPx1 = plab.wf.svd2px(wf_U,sdV);

% (2)
reconstructedFrames = plab.wf.svd2px(wf_U,wf_V); 
sdPx2 = std(reconstructedFrames(:,:,framesToAverageSD),0,3);

% Then plot and use isequal

figure; tiledlayout(1,2);
nexttile;
imagesc(sdPx1);
title('(1) SD of V, then reconstruct pixels');
nexttile;
imagesc(sdPx2);
title('(2) Reconstruct pixels, then take SD');
% Very distinguishable 

isequal(sdPx1, sdPx2); 
dsd = mean(abs(sdPx1 - sdPx2), "all");
fprintf('Mean SD difference is %.12f. \n', dsd);
% Real difference 


% 3) The Bonsai workflow in the loaded data is 'lcr_passive', which is
% passive stimulus presentations as described in the demo 'plab_data_demo'.
% Using what you learned from that demo, make a stimulus-triggered
% (aligned) average movie of widefield fluorescence for -35:35 frames
% around stimulus presentations. Do this separately for stimulus X
% positions = -90 (left), 0 (center) and 90 (right). You'll need to use the
% variable 'stimOn_times', which give the onset time of all stimuli, and
% 'wf_t', which gives the timestamp for each widefield frame.

% So let's build a matrix of X x Y pixel values x +-35f around stimulus
% onsets, with an index for the third dimension of stimulus orientation.
% Then we can use ap. groupfun to make an average movie (I think) and
% ap.imscroll to display it for each one
% actually no, let's not index. the question just wants separately-made
% movies, and that's easier to do

stimOnIdx = interp1(wf_t, 1:numel(wf_t), stimOn_times, 'previous', NaN);

stimOnRange = stimOnIdx + (-35:35);

stimPositionsVec = vertcat(trial_events.values.TrialStimX); % Grab the cells of all stimulus positions

leftRange = stimOnRange(stimPositionsVec == -90,:);
centerRange = stimOnRange(stimPositionsVec == 0,:);
rightRange = stimOnRange(stimPositionsVec == 90,:);

[~, nWindow] = size(leftRange); 

avgLeft = zeros([450,426,nWindow], "double"); % Initialising. 
for i = 1:nWindow
    frameIds = [leftRange(:,i)]; % Get the frame indices for column in question
    leftAvgV = mean(wf_V(:,frameIds),2); % average V selected frames for all components
    leftAvgPx = plab.wf.svd2px(wf_U,leftAvgV);
    avgLeft(:,:,i) = leftAvgPx;
end

avgCenter = zeros([450,426,nWindow], "double"); % Initialising. 
for i = 1:nWindow
    frameIds = [centerRange(:,i)]; % Get the frame indices for column in question
    centerAvgV = mean(wf_V(:,frameIds),2); % average V selected frames for all components
    centerAvgPx = plab.wf.svd2px(wf_U,centerAvgV);
    avgCenter(:,:,i) = centerAvgPx;
end

avgRight = zeros([450,426,nWindow], "double"); % Initialising. 
for i = 1:nWindow
    frameIds = [rightRange(:,i)]; % Get the frame indices for column in question
    rightAvgV = mean(wf_V(:,frameIds),2); % average V selected frames for all components
    rightAvgPx = plab.wf.svd2px(wf_U,rightAvgV);
    avgRight(:,:,i) = rightAvgPx;
end

ap.imscroll(avgLeft);

ap.imscroll(avgCenter);

ap.imscroll(avgRight);

% All looks correct

% I don't really want to do these loops separately but that is what the
% question asked for. Maybe a way to do this with matrices rather than a
% loop? 

% [nComp, ~] = size(wf_V);
% [nTrials, nWindow] = size(leftRange);
% Vcols = wf_V(:, leftRange(:));          
% 
% V3 = reshape(Vcols, nComp, nTrials, nWindow);
% 
% leftAvgV_all = squeeze(mean(V3, 2));   
% 
% leftPxAll = wf_U * leftAvgV_all;       
% 
% avgLeft = reshape(leftPxAll, [450, 426, nWindow]);



% Using "frames" as units to work with widefield data (as in 3 in the
% exercise above, grabbing frames relative to an event) has some downside.
% First, it's less intuitive to use than time, Second, it can be variable
% with respect to an actual event time, because the framerate is
% independent of events - as a diagram, if "|" is a frame capture, "x" is
% an event, and "." is an arbitrary unit of time, this is a possible
% scenario: 
% .....|....x|.....|.....|x....|.....
% Note that the first "x" has a frame in the next timepoint, while the
% second "x" has a frame at a much later relative timepoint. If you just
% look for the next frame after each event "x" then, you're actually
% averaging data from different relative timepoints, slightly smearing your
% data in time. We can account for this by interpolating data between
% frames. For example, in the diagram below, we only have data at "|", but
% we can estimate data at "o" with a linear change in time between frames.
% .....|..o..|.....
% In Matlab, we can do interpolation with the 'interp1' function. This
% function takes a set of known data, and interpolates data between the
% points at specific points. Here's an example, which plots the "actual"
% data, the line between data points, and the interpolated data as red
% points:
% (this interpolates data)
x = [0;1]; % x-values for data
y = [5;6]; % y-values for data
x_interp = [0.2;0.5]; % x-values we want to interpolate y values for
y_interp = interp1(x,y,x_interp); % perform the interpolation
% (this plots the above interpolation)
figure; hold on;
plot(x,y,'.','MarkerSize',20); % plot data
line(x,y); % draw a line for data
plot(x_interp,y_interp,'.r','MarkerSize',20); % plot interpolated data
xlim([-1,2]);
ylim([4,7]);
legend({'Data','Line between data points','Interpolated data'});
xlabel('X');ylabel('Y');

% If our data is multi-dimensional at given timepoints, we can also
% interpolate for each dimension simultanously. Here's an example, note
% it's almost the same as above, except the y-values have 2 dimensions.
x = [0;1]; % x-values for data
y = [5,7;6,8]; % y-values for data
x_interp = [0.2;0.5]; % x-values we want to interpolate y values for
y_interp = interp1(x,y,x_interp); % perform the interpolation
% (this plots the above interpolation)
figure; hold on;
plot(x,y,'.','MarkerSize',20); % plot data
line(x,y); % draw a line for data
plot(x_interp,y_interp,'.r','MarkerSize',20); % plot interpolated data
xlim([-1,2]);
ylim([4,9]);
legend({'Data dim 1','Data dim 2','','','Interpolated data'});
xlabel('X');ylabel('Y');

% Because linear interpolation is a linear operation, for widefield data,
% we can do this on the V's (rather than needing to work with pixels). For
% example, if we want to get a widefield frame at an arbitrary timepoint,
% we can do this:
grab_time = 100; % timepoint to grab (seconds, in Timelite clock)
% (Note: some transposing is required in the interpolation step below: the
% V's are stored as component x time, but 'interp1' requires 'y' data to be
% time x dimension. So, we transpose the input (with '), then transpose the
% output, which make the input the right orientation for 'interp1', and the
% output our normal orientation for V):
V_interp = interp1(wf_t,wf_V',grab_time)';
px_interp = plab.wf.svd2px(wf_U,V_interp);
figure; colormap(gray);
imagesc(px_interp);
axis image
title(sprintf('Frame interpolated at t = %g',grab_time));

% [EXERCISE] 
% Use interpolation to make a stimulus-triggered average as above, but
% using time instead of frames. Use stimulus X position = -90, and grab
% time spanning -1s:1s around the stimulus. Make two movies, one with 0.3s
% resolution (i.e. t = [-1s,-0.7s,-0.4s...]) and one with 0.1s resolution
% (i.e. t = [-1s,-0.9s,-0.8s...])

% An approach that I couldn't quite get to work
% stimTimeRange03 = stimOn_times +(-1:0.3:1);
% VInterp03 = interp1(wf_t,wf_V',stimTimeRange03);
% VInterp03R = reshape(VInterp03,[717,150,7]);
% avgVInterp03R = mean(VInterp03R,2);
% PxInterp03 = squeeze(plab.wf.svd2px(wf_U,avgVInterp03R));
% ap.imscroll(PxInterp03);

% Using a loop instead 

% stimTimeRange03 = stimOn_times +(-1:0.3:1);

stimLeftTimes = stimOn_times(stimPositionsVec ==-90);

res03 = stimLeftTimes +(-1:0.3:1);
res01 = stimLeftTimes +(-1:0.1:1);

nStims = numel(stimLeftTimes);

TLeft03 = zeros([717,numel(res03(1,:)),nStims], "double");
TLeft01 = zeros([717,numel(res01),nStims], "double"); % Initialising. 
for i = 1:nStims
    TLeft03(:,:,i) = interp1(wf_t,wf_V',res03(i))';
end

% --- Operations on widefield SVD spatial components (U's)

% In the same way that we can do linear operations in time on the V's
% above, we can do linear operations in space on the U's. For example, if
% we want to average pixels together in a region-of-interest (ROI), we can
% do that first on the U's and then reconstruct a single trace. 

% Use this code to draw an ROI on the average image:
figure; colormap(gray);
imagesc(wf_avg); 
axis image off
roi_poly = drawpolygon;
% Then use this code to turn your ROI polygon into a binary mask: 
roi_mask = createMask(roi_poly);
figure; imagesc(roi_mask); axis image;
title('ROI mask');

% [EXERCISE] 
% 1) Using the mask 'roi_mask' created above, average the pixel values
% within in the ROI for each component in U, giving you one number for each
% component (rather than a spatial component Y pixels x X pixels x N
% components, you'll have 1 pixel-average x N components). Reconstruct the
% pixel trace with this new averaged spatial component in the same way the
% full-image reconstructions were done above.
%
% 2) Show that this the above method is equivalent to first reconstructing,
% then averaging pixels.

% I have a function 'ap.wf_roi' to make ROIs to make that process easier,
% here's an example use, which can output both the trace and the mask:
[roi_trace,roi_mask2] = ap.wf_roi(wf_U,wf_V,wf_avg);
figure;plot(roi_trace);
title('ROI fluorescence');

% You can also enter a mask as a 5th argument, for example this will use
% exactly the mask you created above, and the trace should be the same as
% the one you reconstructed above:
roi_trace = ap.wf_roi(wf_U,wf_V,wf_avg,[],roi_mask);
figure;plot(roi_trace);
title('ROI fluorescence (with previous mask)');

%% Widefield preprocessing

% The widefield data demoed above has had a few preprocessing steps beyond
% SVD decomposition. Those will be explained here, in the reverse order in
% which they're applied (to work backwards to raw data).

% --- Deconvolution

% Whenever a neuron has a spike (or more readily seen, a burst), calcium
% indicators have a fast rise and then a slow decay. The rise corresponds
% to the spike (note: not the peak, that's delayed from the spike), and the
% slow decay corresponds to nothing - it's an artifact of the indicator.
% This means that the raw fluorescence we see isn't really what we want: we
% only care about the rise. (This is somewhat simplified, for example: if a
% cell is constantly spiking and then stops, this will manifest as a
% baseline fluorescence level which dips, so we care about that drop). The
% fluorescent calcium indicator we primarily use is GCaMP6s. The "s" stands
% for "slow", meaning it's even on the slower end of whats available (but
% much brighter).

% One non-fancy way to approximately get rid of these indicator artifacts
% is to use the derivative of raw fluorescence - this will give us positive
% values on fast rises, and small negative values (which we can ignore) on
% the slow drops. A fancier way is with deconvolution. If we assume that
% each spike produces a characteristic fluorescence rise/fall shape, then
% the fluorescence trace is a convolution of the spike trace with that
% fluorescence shape. We can then deconvolve this rise/fall shape out to
% produce the original spike trace.

% Here's a toy example of the relationship between spikes and fluorescence
% that illustrates the convolution, run this code and zoom into the
% right-hand plot to see the details:
spike_trace = zeros(1000,1); % make an empty trace of spikes
spike_trace(randi(1000,800,1)) = 1; % randomly place spikes
spike_fluorescence = [linspace(0,1,3),linspace(1,0,20)]; % define the fluorescence shape for each spike
fluorescence_trace = conv(spike_trace,spike_fluorescence); % convolve the spike trace with the fluorescence shape
figure;
subplot(1,4,1);
plot(spike_fluorescence);
title('Fluorescence shape for each spike');
subplot(1,4,2:4); hold on;
plot(spike_trace);
plot(fluorescence_trace);
legend({'Spikes','Fluorescence'});

% [EXERCISE] 
% It was mentioned above that taking the derivative is a decent
% approximation for deconvolution. From the above example, plot
% 'spike_trace' with the derivative of 'fluorescence_trace' to see how they
% compare to each other. 

% Deconvolution requires knowing the characteristic fluorescence shape (aka
% "kernel") for each spike. We derived this data from simultaneous imaging
% and electrophysiology, and finding the best kernel to translate
% fluorescence into multiunit spiking (Peters et al. 2021, Extended Data
% Fig. 4a-b: https://www.nature.com/articles/s41586-020-03166-8/figures/9).

% We apply that deconvolution kernel with the function 'ap.wf_deconv'.
% Because deconvolution is a linear operation, it can be done directly on
% the V's from SVD. In the loaded example dataset, non-deconvolved data is
% stored as 'wf_Vdf';

% [EXERCISE] 
% 1) Reconstruct pixels from non-deconvolved ('wf_Vdf') and deconvolved
% ('wf_V') from frames 150-300. Normalize the values in each movie (e.g.
% make the range of values from 0-1, because they have very different
% scales), horizontally concatenate the movies, and scroll through them
% with 'ap.imscroll'. What features are different between non-deconvolved
% and deconvolved fluorescence?
%
% 2) In the previous section, you made stimulus-triggered average widefield
% movies. Do that again with with stimulus X = -90 using non-deconvolved
% data, and compare to deconvolved data. There should be a region in the
% visual cortex with responds to the stimulus - draw an ROI over this
% region for non-deconvolved and deconvolved data, and plot those traces
% (normalized to have similar ranges) on top of each other. (One easy way to
% draw an ROI on a movie and get the associated trace: view a
% movie with 'ap.imscroll' and press 'r' on the figure, which lets you draw a
% polygon. After completing and double-clicking the polygon, the structure
% 'roi' will be saved in the workspace, which includes fields roi.trace
% (the average values within the ROI) and roi.mask (the binary mask for
% your ROI)). 

% --- Normalization (dF/F)

% Absolute values of fluorescence don't have any intrinsic information. For
% example: if we have our LED dim one day and bright the next, the
% fluorescence brightness will change, but this is unrelated to activity.
% Similarly, if area 1 has more GCaMP than area 2, it will be brighter, but
% this is not because it is more active. Instead, only relative values
% carry information about activity, specifically: how does fluorecence
% change relative to "baseline"? We therefore want to "normalize" the
% fluorescence in each pixel to its baseline.

% The standard way to normalize fluorescence is as the difference between
% measured fluorescence (F) and baseline fluorescence (F_0), divided by
% F_0, = (F-F_0)/F_0, = ΔF/F_0. Our resulting units then are "fractions of
% baseline fluorescence": at baseline = 0, double baseline = 1, half
% baseline = -1.

% What we consider "baseline" can vary - ideally it's when there's "nothing
% happening", but we just use the average fluorescence across the day as a
% reasonable and easy reference. Have a look at the average image in
% pseudocolor to highlight differences, you can see that there's a range of
% brightness across the brain:
figure;
imagesc(wf_avg);
axis image;
colormap(hsv);

% [EXERCISE] 
% Normalized data is loaded as 'wf_Vdf', and non-normalized data is loaded
% as 'V_neuro_hemocorr'. Reconstruct frames 580-600 from non-normalized and
% normalized data and view the movies - there should be a spots of activity
% in the lower left (visual cortex) and center left (motor cortex). Get
% activity in ROIs over those areas for both non-normalized and normalized
% data. Make two plots with overlaid visual/motor traces: one for
% non-normalized data, and one with normalized data. What's the difference
% between these plots? Why is that difference present?

% --- Hemodynamic correction

% Blood in the brain introduces an artifact in our imaging: blood is routed
% to active areas, which obscures the fluorescence and darkens the signal.
% This is generally slow (starts ~0.5s after activity and can last ~2s),
% and is almost removed just by deconvolution (which pulls out only fast
% components), but we still want to remove this artifact to have more
% accurate data.

% The technique to remove this artifact is to use two excitation lights,
% blue and violet, to trigger light emission from GCaMP in different ways.
% When GCaMP is excited by blue light, it emits more photons if it's bound
% to calcium than if it's not (which is why we get activity-dependent
% fluorescence). When GCaMP is excited by violet light, it emits the same
% amount of photos regardless of whether it's bound to calcium or not (so
% the violet-excited GCaMP signal is activity-independent).

% We then have activity-dependent (blue) and activity-independent (violet)
% signals, both of which are sometimes dimmed by blood. If we properly
% scale the violet signal and subtract it from the blue signal, we then
% remove the blood (hemodynamic) effect, leaving just the activity effect.

% This is corrected in the loading script, but we can still look at these
% components separately. So far, you've been using the
% hemodynamically-corrected data (Udf, fVdf), where the uncorrected
% versions are in:
% blue (activity): Un, fVn (n for 'neural')
% violet (hemodynamic): Uh, fVh (h for 'hemodynamic');

% [EXERCISE] 
% Non-hemodynamic-corrected (raw) data is loaded as 'wf_V_raw{1}'. As
% before, make a stimulus-triggered widefield movie for stimulus X = -90
% spanning -0.5s:+4s, for both raw data and hemodynamic-corrected data
% ('V_neuro_hemocorr'). Draw an ROI over the responsive areas in both
% movies and plot the traces for both ROIs together on one plot. What are
% the differences between the traces?

%% Widefield responses across the cortex

% Load this data for the following section: 
animal = 'AP010';
rec_day = '2023-08-13';
rec_time = '1555';
verbose = true;
ap.load_recording;

% --- Regional correlations in widefield data

% Our widefield images span the whole dorsal sorface of the brain. To get a
% feel for where areas are, it can be useful to look at correlations
% between pixels. This is a tool to explore that: it gives a map of
% correlations between a clicked pixel and all other pixels. Try clicking
% around the brain to get a sense for the pattern of correlations, also you
% can press 'h' to enable 'hover mode' that follows your cursor:
ap.wf_corrviewer(wf_U,wf_V);

% [EXERCISE] 
% In the correlation viewer above: 
%
% 1) There should be some regions that have topographic relationships,
% where adjacent pixels in one region are correlated to adjacent pixels in
% another region. You should see the hotspots of correlation move in an
% ordered way across multiple regions, and then come together at specific
% points. Which regions have this relationship? Where do the correlations
% join together? 
% 
% 2) If you had to divide the cortex into regions based on the correlation
% patterns, where would you draw the lines? Take a look at Figure 3 in this
% paper which labels dorsal cortical areas:
% https://www.cell.com/cell/fulltext/S0092-8674(20)30402-5.
% What regions do your divisions correspond to?

% --- Regional activity during events

% In the loaded dataset, the mouse is performing a task where a stimulus
% is presented on the right-hand screen, then moving the wheel brings the
% stimulus to the center to trigger a sucrose reward. 

% This means there are three main things that happen in this task: 
% 1) the stimulus comes on (with the time stored in 'stimOn_times')
% 2) the movement starts (the first movement after the stimulus can be
% found using the 'wheel_move' variable mentioned in 'plab_data_demo')
% 3) the reward is delivered and consumed (reward delivery times are stored
% in the variable 'reward_times')

% [EXERCISE] 
% 1) Make three event-triggered averages of widefield activity for each of
% the 3 events above: stimulus onset, onset of movement after the stimulus
% is presented, reward delivery. Comparing these movies to the mouse atlas,
% which region are active during these events?
%
% 2) The "reaction time" (defined here as the time from stimulus onset to
% movement) varies on each trial. Plot the activity in the visual cortex
% and motor cortex on each trial, sorted by reaction time, by following
% these steps. First, get the activity in an ROI over the visual and motor
% cortex. Next, align data from each ROI to stimulus onsets. Calculate the
% reaction time for each trial, then plot a heatmap (using 'imagesc') of
% the data from each ROI, with time on the x-axis and trial sorted by
% reaction time on the y-axis. Draw a line at t = 0 for reference, and
% another at the movement onset time (note this will be different for each
% trial). What do these plots show about alignment of activity to these
% events?

% --- Retinotopy of visual areas

% Visual regions of the cortex can be mapped by showing small squares
% across the screen in a sparse, uncorrelated pattern ("sparse noise"), and
% finding which places on the cortex respond to squares at each position.

% Each visual region has a characteristic "visual field sign", which is the
% orientation of the visual representation with response to the world: a
% positive sign means the orientation is the same as the world (up on the
% screen is anterior in the brain, right on the screen is right on the
% brain), a negative sign means the orientation is the mirror of the world
% (up on the screen is posterior on the brain, right on the screen is left
% on the brain)

% We can also use visual field sign to align our widefield images to the
% mouse atlas (Allen Common Coordinate Framework, or "CCF"), because we
% know where the regions should be and what sign they should have.

% The function 'plab.wf.retinotopy_vfs' loads data for a recording and
% calculates the visual field sign (vfs). Run the following command to get
% and plot the vfs for an example recording: 
animal = 'AP010';
rec_day = '2023-08-09';
rec_time = '1736';
vfs = plab.wf.retinotopy_vfs(animal,rec_day,rec_time);

% [EXERCISE]
% Compare the regions from the visual field sign to the pixel correlation
% using 'ap.wf_corrviewer'. Try to identify visual regions on the
% correlation viewer using patterns of correlation and referring to the
% visual field sign map - how does the visual field sign show up in the
% pattern of correlations?

%% Next demo

% To learn how to align data across days and animals (e.g. to look for
% changes across days, to average across animals), move on to
% plab_widefield_demo2







