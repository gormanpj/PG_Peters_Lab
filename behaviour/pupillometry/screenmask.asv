% Code to trim out mousecam frames taken when the stimscreen is on 

% Once recording is loaded

mousecamOnIdx = find(diff([double(mousecam_thresh)]) == 1) + 1 ;  
mousecamOffIdx = find(diff([double(mousecam_thresh)]) == -1) + 1 ;  

screenOnFrames = [];
screenOnProportion = [];

% Exposure data contains instances that don't make it into final frames.
% This makes a mask to get rid of those 
finalExposesMask = mousecam_exposeOn_times >= mousecam_times(1) & mousecam_exposeOn_times <= mousecam_times(end);

for n = 1:numel(mousecamOffIdx);
    exposeMask = (mousecamOnIdx(n):1:(mousecamOffIdx(n) -1));
    screenOnFrames(n) = sum(screen_on(exposeMask));
end

screenOff = screenOnFrames' == 0;
screenOffMask = screenOff(finalExposesMask);

figure; histogram(screenOnFrames);