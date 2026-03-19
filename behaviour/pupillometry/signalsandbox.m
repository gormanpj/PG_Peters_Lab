sgbase = sgolayfilt(diameterPx_all{4} ,3,9);
figure; hold on
plot(diameterPx_all{4}, 'k')
plot(sgbase, 'r')







base = diameterPx_all{4};
basefill = fillmissing(base, "linear");
orig_mean = mean(basefill);
figure;hold on 
plot(basefill,'k');
fRanges = [9 9.5; 11.5 11.75]; % 2.4 2.5; 6.73 6.8; 

for i = 1:size(fRanges, 1)
    basefill = bandstop(basefill, fRanges(i,:), 30);
end
dc_restored_basefill = basefill + (orig_mean - mean(basefill));

plot(basefill, 'c');


f = fRanges(1,:);                 % test first stopband only
[y,d] = bandstop(basefill, f, Fs);% returns digitalFilter d
disp(d)
fvtool(d)


% 1) basic stats
orig_mean = mean(basefill);
y = basefill; for i=1:size(fRanges,1); y = bandstop(y, fRanges(i,:), Fs); end
filt_mean = mean(y);
fprintf('mean(orig)=%.6g, mean(filt)=%.6g, diff=%.6g\n', orig_mean, filt_mean, filt_mean-orig_mean);

% 2) inspect the filter (single-band test)
f = fRanges(1,:);
[y1,d] = bandstop(basefill, f, Fs);  % get digitalFilter object d
try
    % try to get transfer function coefficients if possible
    [b,a] = tf(d);       % works for many digitalFilter objects
catch
    % fallback: visualize in fvtool
    disp('Cannot extract b,a with tf(d). Opening fvtool(d).')
    fvtool(d); return
end

% 3) DC gain from b,a
dcGain = sum(b)./sum(a);   % DC gain = H(z=1)
fprintf('DC gain for test filter = %.6g\n', dcGain);

% 4) plot the transfer function around 0..Fs/2
freqz(b,a,2048,Fs);
title('Filter magnitude / phase');

% confirm Fs and signal length
fprintf('Fs=%g, length=%d\n', Fs, length(basefill));

% PSD (better windowing/resolution)
figure; pwelch(basefill, hann(2048), 1024, 8192, Fs);
title('PSD (pwelch)');

% Spectrogram (time-varying view)
figure;
window = 1024; noverlap = 768; nfft = 2048;
[s,f,t,p] = spectrogram(basefill, window, noverlap, nfft, Fs, 'yaxis');
imagesc(t,f,10*log10(abs(p))); axis xy;
xlabel('Time (s)'); ylabel('Freq (Hz)'); colorbar; title('Spectrogram (dB)');


Fs = 30;         % confirm
y = basefill;

% 1) Notch each spike with moderate Q
for i = 1:size(fRanges,1)
    f0 = mean(fRanges(i,:));
    bw = fRanges(i,2) - fRanges(i,1);
    if bw <= 0, continue; end
    Q = max(15, f0/bw);            % choose Q; cap to avoid extreme Q
    w0 = f0/(Fs/2);
    if w0 <= 0 || w0 >= 1, continue; end
    [bn,an] = iirnotch(w0, w0/Q);
    y = filtfilt(bn,an,y);
end

% 2) Savitzky-Golay smoothing to remove residual high-frequency noise
% choose odd frame length in samples. e.g. frame=31, polyorder=3
frame = 31; polyorder = 3;
frame = min(frame, length(y)-1 + mod(length(y)-1,1)); % sanity-check
y_sg = sgolayfilt(y, polyorder, frame);

figure; hold on;
plot(basefill, 'k'); plot(y, 'c'); plot(y_sg, 'r'); legend('orig','notched','notched+SG');



% design lowpass (zero-phase)
Fc = 5;              % choose cutoff in Hz (set to suit your signal)
Wn = Fc/(Fs/2);
[b,a] = butter(4, Wn, 'low');
y_lp = filtfilt(b,a,basefill);

% then Savitzky-Golay for small residual smoothing
y_lp_sg = sgolayfilt(y_lp, 3, 31);

figure; hold on;
plot(basefill,'k'); plot(y_lp,'c'); plot(y_lp_sg,'r'); legend('orig','lowpass','lowpass+SG');



% lowpass then SG
base = diameterPx_all{4};
basefill = fillmissing(base, "linear");
orig_mean = mean(basefill);
figure;hold on 
plot(basefill,'k');

%lowpass is causal
baselow = lowpass(fillmissing(basefill, "linear"), 4, 30);
plot(baselow, 'g')

% SG is not
baselowsg = sgolayfilt(baselow, 3, 15);
plot(baselowsg, 'y')



% lowpass then movmean
base = diameterPx_all{4};
basefill = fillmissing(base, "linear");
orig_mean = mean(basefill);
figure;hold on 
plot(basefill,'k');

%lowpass is causal
baselow = lowpass(fillmissing(basefill, "linear"), 4, 30);
plot(baselow, 'g')

%then a causal movmean
baselowmean = movmean(baselow, [3 0]);
plot(baselowmean, 'b')








%%% Blues and violets

violetOn = timelite.data(:,8) >= 0.05;
blueOn = timelite.data(:,9) >= 0.05;

for k = 1:numel(mousecam_times)
    exposeMask = isbetween(timelite.timestamps, mousecam_times(k), mousecamOffTimes(k));
    violetOnFrames(k) = any(violetOn(exposeMask));
    violetOnFramesCt(k) = sum(violetOn(exposeMask));
    end
figure;histogram(violetOnFramesCt);
for k = 1:numel(mousecam_times)
    exposeMask = isbetween(timelite.timestamps, mousecam_times(k), mousecamOffTimes(k));
    blueOnFrames(k) = any(blueOn(exposeMask));
    blueOnFramesCt(k) = sum(blueOn(exposeMask));
end
figure;histogram(blueOnFramesCt);

violetOnMask = violetOnFramesCt >= 3;
blueOnMask = blueOnFramesCt >= 3;

blues = diameterPx_all{4};
violets = diameterPx_all{4};
blues(~blueOnMask) = NaN;
violets(~violetOnMask) = NaN;
figure; hold on
plot(diameterPx_all{4},'k')
plot(blues,'b')
plot(violets,'m')

bluesfilled = fillmissing(blues, "linear");
violetsfilled = fillmissing(violets, "linear");
plot(diameterPx_all{4},'k')
plot(bluesfilled,'b')
plot(violetsfilled,'m')

indigos = (bluesfilled+violetsfilled)/2;
plot(indigos, 'c')

indigoLow = lowpass(indigos, 6, 30);
plot(indigoLow, 'r');

rights = psthData(orientationIdx == 1,:);

figure; hold on 
for stims = 1:50
    plot(rights(stims,:), 'k')
end

avg = mean(rights, 1,"omitmissing");
plot(avg, 'r', 'LineWidth',2)

figure; tiledlayout
for circ = 1:nRecs
    rights = psthData(orientationIdx == 3 & dayIdx == circ,:);
    
    nexttile; hold on 
    for stims = 1:50
        plot(rights(stims,:), 'k')
    end
    avg = mean(rights, 1,"omitmissing");
    plot(avg, 'r', 'LineWidth',2)
    hold off
end


figure;hold on 
for toodle = 1:9
    plot(diameterPx_all{toodle})
end 
hold off

rights = psthData_all(orientationIdx_all == 3 & learnDayIdx_all == 1,:);

figure; tiledlayout 
for day = min(learnDayIdx_all):max(learnDayIdx_all)
    rights = psthData_all(orientationIdx_all == 3 & learnDayIdx_all == day,:);
    nexttile; hold on 
    for stims = 1:size(rights,1)
        plot(rights(stims,:), 'k')
    end
    avg = mean(rights, 1,"omitmissing");
    plot(avg, 'r', 'LineWidth',2)
    title(day)
    xline(21, 'k--')
    hold off
end

col = ['b', 'k', 'r'];
for day = -2:1 
    figure; hold on 
    for orientation = 1:3
        rights = psthData_allQui(orientationIdx_allQui == orientation & learnDayIdx_allQui == day,:);
        avg = mean(rights, 1,"omitmissing");
        bline = avg(21);
        avgsub = avg-bline;
        plot(avgsub, col(orientation), 'LineWidth',2)
    end
    xline(21, 'k--')
    title(day)
    hold off
end



orientNames = {'left','center','right'};

figure;
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

for oi = 1
    figure; hold on
    
    % All trials for this orientation
    allMask = orientationIdx_all == oi & learnDayIdx_all == -1;
    allTraces = psthData_all(allMask, :);
    
    % Quiescent trials for this orientation
    quiMask = orientationIdx_allQui == oi & learnDayIdx_allQui == -1;
    quiTraces = psthData_allQui(quiMask, :);
    
    % Plot all trials in gray
    for r = 1:size(allTraces,1)
        plot(x, allTraces(r,:), 'Color', [0.75 0.75 0.75]);
    end
    
    % Mean lines
    if ~isempty(allTraces)
        plot(x, mean(allTraces, 1, 'omitnan'), 'k', 'LineWidth', 2);
    end
    if ~isempty(quiTraces)
        plot(x, mean(quiTraces, 1, 'omitnan'), 'b', 'LineWidth', 2);
    end
    
    title(orientNames{oi});
    xlabel('Frame');
    ylabel('Diameter');
    %legend({'all trials','all mean','quiescent mean'}, 'Location', 'best');
    hold off
end

figure; hold on

% Plot all trials in light gray
for r = 1:size(psthData_all,1)
    plot(x, psthData_all(r,:), 'Color', [0.8 0.8 0.8]);
end

% Mean of all trials
plot(x, mean(psthData_all, 1, 'omitnan'), 'k', 'LineWidth', 2);

% Mean of quiescent trials
plot(x, mean(psthData_allQui, 1, 'omitnan'), 'g', 'LineWidth', 2);

xlabel('Frame');
ylabel('Diameter');
title('All trials');
%legend({'all trials','all mean','quiescent mean'}, 'Location', 'best');

hold off