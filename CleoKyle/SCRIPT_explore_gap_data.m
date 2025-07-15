%%
addpath(genpath('C:\Users\mattm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath('C:\Users\mattm\Documents\GitHub\nsb-mouse-2025\CleoKyle');

cfg_master = [];
cfg_master.debug = 0;
cfg_master.frame_pulse_thr = 1.5;
cfg_master.ybounds = [400 650]; % restrict video samples to this range
cfg_master.xbounds = [500 1750];
cfg_master.gap = 'landing'; % 'takeoff', 'midpoint', 'landing'; uses edges.csv file
cfg_master.norm = 'zscore'; % 'dF', 'zscore'
cfg_master.trackPoint = 1; % 2 is mid-body, 1 is nose
cfg_master.writeFig = 1; % write figure image files?
cfg_master.peth_normalize = []; % [] (no normalization), 'range', 'zscore'

%% load
fname = FindFile('*.h5'); % SLEAP analysis output
occupancy_matrix = h5read(fname,'/track_occupancy');
tracks_matrix = h5read(fname,'/tracks');

video_n_frames = length(occupancy_matrix);
fprintf('SLEAP output has %d frames\n', video_n_frames);

%% SLEAP output is indexed by video 1 frames, not seconds, so need to convert
fname = FindFile('acq*.csv'); % Digilent raw data file. At NS&B 2025, Ch1 is GRAB-DA, Ch2 is video 1 frame trigger

DATA = readtable(fname);
DATA = DATA(2:end,:); % first sample can be artifact, so just remove

DATA.Time_s_ = DATA.Time_s_ - DATA.Time_s_(1);

if cfg_master.debug
    figure(1)
    ax1 = subplot(311); plot(DATA.Time_s_, DATA.Channel1_V_); title('Raw fiber signal')
    ax2 = subplot(312); plot(DATA.Time_s_,DATA.Channel2_V_); title('Ch2 Video 1 camera pulse')
end

Ch2_diff = diff(DATA.Channel2_V_ > cfg_master.frame_pulse_thr);
Ch2_up_ix = find(Ch2_diff > 0); Ch2_frame_time_sec = DATA.Time_s_(Ch2_up_ix);

video_tvec = DATA.Time_s_(Ch2_up_ix); video_n_pulses = length(Ch2_frame_time_sec);
fprintf('Digilent data has %d pulses\n', video_n_pulses);

if video_n_pulses > video_n_frames
    video_tvec = video_tvec(1:video_n_frames); % this gives time for every frame
end
video_x = tracks_matrix(:,cfg_master.trackPoint,1);
video_y = tracks_matrix(:,cfg_master.trackPoint,2);

if cfg_master.debug, subplot(313), plot(video_x, video_y, '.k'); axis off; hold on; end

%% load the trial data
ReadObstacleExpKeys; nT = length(start_frames);
trial_start = video_tvec(start_frames); trial_end = video_tvec(end_frames);

edges = readtable('edges.csv');

cfg_master.gapx = [];
switch cfg_master.gap % load edges and find midpoint
    case 'midpoint' % mean of edges
    for iT = 1:nT
        this_row = find(edges.Var1 == iT);
        cfg_master.gapx(iT) = mean([edges(this_row,:).Var2 edges(this_row,:).Var3]);
    end
    case 'takeoff'
        for iT = 1:nT
            this_row = find(edges.Var1 == iT);
            cfg_master.gapx(iT) = edges(this_row,:).Var2;
        end
    case 'landing'
        for iT = 1:nT
            this_row = find(edges.Var1 == iT);
            cfg_master.gapx(iT) = edges(this_row,:).Var3;
        end
    otherwise
        error('Unknown cfg_master.gapx')
end

%% preprocess the fiber data
Fs = 1 ./ median(diff(DATA.Time_s_));  % Sampling frequency in Hz
fprintf('Fs = %.2f Hz\n', Fs)

data = detrend(DATA.Channel1_V_, 'linear', 2); % whole session detrending

Fc = 20;  % Cutoff frequency in Hz
[b, a] = butter(4, Fc / (Fs / 2), 'low');
data = filtfilt(b, a, data);

switch cfg_master.norm
    case 'zscore'
        fiber_data = (data - mean(data)) / std(data); % z-score
        fiber_tvec = DATA.Time_s_;
    case 'dF'
        % moving window
        winsize = 5000;
        nS = length(data);
        data_dF = nan(nS,1);
        fprintf('Computing dF/F...\n');
        for iS = 1:nS

            start = max(1, iS - winsize - 1);
            this_window_data = data(start:iS-1); % window up until this data point
            data_dF(iS) = (data(iS) - mean(this_window_data)) ./ std(this_window_data);

        end
        fiber_data = data_dF;
        fiber_tvec = DATA.Time_s_;
    otherwise
        error('Unknown cfg_master.norm method')
end

%% preprocess the video data & find gap crossings
video_keep = (video_x >= cfg_master.xbounds(1) & video_x <= cfg_master.xbounds(2) ...
    & video_y >= cfg_master.ybounds(1) & video_y <= cfg_master.ybounds(2));
video_x(~video_keep) = NaN; video_y(~video_keep) = NaN;

if cfg_master.debug, figure(1); subplot(313); plot(video_x, video_y, 'or'); end

nT = length(trial_start);
gap_times = nan(nT,1);
for iT = 1:nT

    this_trial_idx = find(video_tvec >= trial_start(iT) & video_tvec < trial_end(iT));
    this_trial_tvec = video_tvec(this_trial_idx);
    this_trial_x = video_x(this_trial_idx); this_trial_y = video_y(this_trial_idx);
    
    if cfg_master.debug, figure(2)
        subplot(211); plot(this_trial_x, this_trial_y, 'b'); hold on;
        title(sprintf('Position filtering check, trial %d',iT));
        subplot(212); plot(this_trial_tvec, this_trial_x), 'b'; hold on;
    end

    this_trial_x = medfilt1(this_trial_x, 9, 'omitnan', 'truncate'); this_trial_y = medfilt1(this_trial_y, 9, 'omitnan', 'truncate');
    
    if cfg_master.debug
        subplot(211); plot(this_trial_x, this_trial_y, 'r'); hold on;
        subplot(212); plot(this_trial_tvec, this_trial_x, 'r'); hold on;
    end

    good_idx = find(~isnan(this_trial_x) & ~isnan(this_trial_y));
    x_interp = interp1(this_trial_tvec(good_idx), this_trial_x(good_idx), this_trial_tvec, 'linear');
    y_interp = interp1(this_trial_tvec(good_idx), this_trial_y(good_idx), this_trial_tvec, 'linear');

    if cfg_master.debug
        subplot(211); plot(x_interp, y_interp ,'g'); hold off;
        subplot(212); plot(this_trial_tvec, x_interp, 'g'); hold off;
        pause;
    end

    video_x(this_trial_idx) = x_interp; 
    video_y(this_trial_idx) = y_interp;

    % find gap crossing if it exists
    crossed_gap = x_interp < cfg_master.gapx(iT);
    crossed_gap = medfilt1(double(crossed_gap), 9, 'omitnan', 'truncate');
    temp = diff(crossed_gap);
    cross_idx = find(temp == 1);

    switch length(cross_idx)
        case 0
            fprintf('No crossing found on trial %d.\n',iT)
        case 1
            fprintf('Gap crossing found on trial %d.\n', iT);
            gap_times(iT) = this_trial_tvec(cross_idx);
        otherwise
            fprintf('Multiple crossings (n = %d) for trial %d, using only first one...\n', length(cross_idx), iT);
            gap_times(iT) = this_trial_tvec(cross_idx(1));
    end


end

% get velocity & acceleration
xdot = dxdt(video_tvec, video_x); ydot = dxdt(video_tvec, video_y);
v = sqrt(xdot.^2 + ydot.^2); vdot = dxdt(video_tvec,v);

%%
figure(2); 
xa = subplot(311);
xline(trial_start,'g--', 'LineWidth',2); hold on; xline(trial_end, 'r--', 'LineWidth',2); xline(gap_times, 'm--', 'LineWidth',2);
yline(cfg_master.gapx,'m:','LineWidth',2);
plot(video_tvec, video_x, '.k', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
ylabel('position (pix)')
[~, fd, ~] = fileparts(pwd); th = title(fd); set(th, 'Interpreter', 'none')

fa = subplot(312);
xline(trial_start,'g--', 'LineWidth',2); hold on; xline(trial_end, 'r--', 'LineWidth',2); xline(gap_times, 'm--', 'LineWidth',2);
plot(fiber_tvec, fiber_data, 'g', 'LineWidth',2);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
ylabel(sprintf('DA (%s)', cfg_master.norm));

vdota = subplot(313);
plot(video_tvec, xdot, 'r', 'LineWidth', 2);
hold on;
xline(trial_start,'g--', 'LineWidth',2); xline(trial_end, 'r--', 'LineWidth',2); xline(gap_times, 'm--', 'LineWidth',2);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
ylabel('velocity')
yl = ylim; ylim([-200 yl(2)])

linkaxes([xa fa vdota], 'x');
maximize

%% loop over trials to output figures
for iT = 1:length(trial_start)

    subplot(311)
    xlim([trial_start(iT)-5 trial_end(iT)+5]);
    title([fd ' trial ' num2str(iT)]);

    drawnow;
    if cfg_master.writeFig
        WriteFig(gcf, ['trial ' num2str(iT) '_'  cfg_master.gap num2str(cfg_master.trackPoint) '_' cfg_master.norm], 'pngOnly');
    end
end
close(gcf)

%% make PETHs
fiber_tsd = tsd(fiber_tvec, fiber_data');
cspace = linspace(0.2,1,nT);

% start-fiber
cfg_peth = []; cfg_peth.normalize = cfg_master.peth_normalize;
[startPETH, startPETHt] = TSDpeth(cfg_peth, fiber_tsd,  trial_start);

subplot(221)
xline(0,'g--', 'LineWidth',1);
for iT = 1:nT, plot(startPETHt.tvec, startPETHt.data(iT,:), 'Color', [0.2 cspace(iT) 0.2], 'LineWidth', 0.5); hold on; end

plot(startPETH, 'g', 'LineWidth', 3);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('DA (%s)', cfg_master.norm));
title('trial start')
xlim([-2 2]);

% start-velocity
vel_tsd = tsd(video_tvec, xdot);

cfg_peth = [];
[startPETH, startPETHt] = TSDpeth(cfg_peth, vel_tsd, trial_start);
subplot(223)
xline(0,'g--', 'LineWidth',1);
for iT = 1:nT, plot(startPETHt.tvec, startPETHt.data(iT,:), 'Color', [cspace(iT) 0.2 0.2], 'LineWidth', 0.5); hold on; end

plot(startPETH, 'r', 'LineWidth', 3);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel('velocity');

% cross-fiber
gap_timesnn = gap_times(~isnan(gap_times)); nGnn = length(gap_timesnn);

cfg_peth = []; cfg_peth.normalize = cfg_master.peth_normalize;
[gapPETH, gapPETHt] = TSDpeth(cfg_peth, fiber_tsd, gap_timesnn);

subplot(222)
xline(0,'m--', 'LineWidth',1);
for iT = 1:nGnn, plot(gapPETHt.tvec, gapPETHt.data(iT,:), 'Color', [0.2 cspace(iT) 0.2], 'LineWidth', 0.5); hold on; end

plot(gapPETH, 'g', 'LineWidth', 3);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('DA (%s)', cfg_master.norm));
title('gap crossing')
xlim([-2 2]);

% cross-velocity
vel_tsd = tsd(video_tvec, xdot);

cfg_peth = [];
[gapPETH, gapPETHt] = TSDpeth(cfg_peth, vel_tsd, gap_timesnn);
subplot(224)
xline(0,'m--', 'LineWidth',1);
for iT = 1:nGnn, plot(gapPETHt.tvec, gapPETHt.data(iT,:), 'Color', [cspace(iT) 0.2 0.2], 'LineWidth', 0.5); hold on; end

plot(gapPETH, 'r', 'LineWidth', 3);
set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel('velocity');

maximize;
%%
if cfg_master.writeFig
    WriteFig(gcf, ['PETHs_'  cfg_master.gap num2str(cfg_master.trackPoint) '_' cfg_master.norm], 'pngOnly');
    close(gcf)
else
    figure
end
%% as above but now by trial type
cols = 'rgb';
tt = unique(trial_types); if length(tt) > 3, tt = tt(1:3); end % allow only 3 trial types
cfg_peth = []; cfg_peth.normalize = cfg_master.peth_normalize;

% start-fiber
subplot(221)
xline(0,'g--', 'LineWidth',1);

clear h
for iTT = 1:length(tt)
    this_idx = find(trial_types == tt(iTT));
    [startPETH, startPETHt] = TSDpeth(cfg_peth, fiber_tsd, trial_start(this_idx));
    h(iTT) = plot(startPETH, 'Color', cols(iTT), 'LineWidth', 3); hold on;
end

set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('DA (%s)', cfg_master.norm));
title('trial start')
xlim([-2 2]);
%legend(h, num2str(tt(:)))

% start-velocity
subplot(223)
xline(0,'g--', 'LineWidth',1);

cfg_peth = [];
clear h
for iTT = 1:length(tt)
    this_idx = find(trial_types == tt(iTT));
    [startPETH, startPETHt] = TSDpeth(cfg_peth, vel_tsd, trial_start(this_idx));
    h(iTT) = plot(startPETH, 'Color', cols(iTT), 'LineWidth', 3); hold on;
end

set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('velocity'));
title('trial start')
xlim([-2 2]);
legend(h, num2str(tt(:)), 'Location', 'northwest'); legend boxoff;

% cross-fiber
subplot(222)
xline(0,'m--', 'LineWidth',1);
cfg_peth = []; cfg_peth.normalize = cfg_master.peth_normalize;

clear h
for iTT = 1:length(tt)
    this_idx = find(trial_types == tt(iTT));
    this_gap = gap_times(this_idx);
    if ~all(isnan(this_gap))
        [startPETH, startPETHt] = TSDpeth(cfg_peth, fiber_tsd, this_gap(~isnan(this_gap)));
        h(iTT) = plot(startPETH, 'Color', cols(iTT), 'LineWidth', 3); hold on;
    end
end

set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('DA (%s)', cfg_master.norm));
title('gap crossing')
xlim([-2 2]);

% cross-velocity
subplot(224)
xline(0,'m--', 'LineWidth',1);
cfg_peth = [];

clear h
for iTT = 1:length(tt)
    this_idx = find(trial_types == tt(iTT));
    this_gap = gap_times(this_idx);
    if ~all(isnan(this_gap))
        [startPETH, startPETHt] = TSDpeth(cfg_peth, vel_tsd, this_gap(~isnan(this_gap)));
        h(iTT) = plot(startPETH, 'Color', cols(iTT), 'LineWidth', 3); hold on;
    end
end

set(gca, 'TickDir', 'out', 'LineWidth', 1, 'FontSize', 18); box off; grid on;
xlabel('time (s)'); ylabel(sprintf('velocity'));
title('gap crossing')
xlim([-2 2]);

maximize;

%%
if cfg_master.writeFig
    WriteFig(gcf, ['PETHs_byTrial_' cfg_master.gap num2str(cfg_master.trackPoint) '_' cfg_master.norm], 'pngOnly');
    close(gcf)
end
