% Synchronize the video with the mice.
% Cowen 2024
% MvdM 2025 adapted for 2 video streams
% analog ch2: video stream 1
% DIO 0: video stream 2

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% some options and params
cfg = [];
cfg.trials_only = 1; % extract trial info from ExpKeys
cfg.trial_padding = 25; % number of frames either side of start and end frames

% define the files.
close all
fd = 'C:\data\nsb2025\M288-2025-07-09_gap';
cd(fd);
digilent_file = [fd '\acq0020.csv'];
video_file{1} = [fd '\video1_qtr.mp4'];
video_file{2} = [fd '\video2_qtr.mp4'];
dlc_file = []; % The deep lab cut csv file if it exists
out_video_dir = "C:\Temp\";
vid_limits_x = [0 960];
vid_limits_y = [0 600];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,digilent_fname] = fileparts(digilent_file);
out_video_file = fullfile(out_video_dir,digilent_fname);

% load the digilent data and get the timestamps for each video frame. 
DATA = readtable(digilent_file);
DATA = DATA(2:end,:); % first sample can be artifact

DATA.Time_s_ = DATA.Time_s_ - DATA.Time_s_(1);

Fs = 1000;  % Sampling frequency in Hz
data = detrend(DATA.Channel1_V_, 'linear'); % whole session detrending

Fc = 20;  % Cutoff frequency in Hz
[b, a] = butter(4, Fc / (Fs / 2), 'low');
data = filtfilt(b, a, data);

DATA.smoothed = (data - mean(data)) / std(data);

figure
ax1 = subplot(311);
plot(DATA.Time_s_, DATA.smoothed);
title('Fiber signal')

ax2 = subplot(312);
plot(DATA.Time_s_,DATA.Channel2_V_)
title('Analog Ch2')

ax3 = subplot(313);
plot(DATA.Time_s_,DATA.DIO0)
title('DIO0')

linkaxes([ax1 ax2 ax3], 'x');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate and compute sample rates and num frames.

Ch2_diff = diff(DATA.Channel2_V_ > 1.5);
Ch2_up_ix = find(Ch2_diff > 0);
Ch2_frame_time_sec = DATA.Time_s_(Ch2_up_ix);
Ch2_n_frames = length(Ch2_frame_time_sec);
    Ch2_frame_interval_sec = median(diff(Ch2_frame_time_sec));
Ch2_video_sFreq = 1/Ch2_frame_interval_sec;

fprintf('Analog Ch2: %d frames, %1.2f Hz\n',Ch2_n_frames, Ch2_video_sFreq)

if (DATA.Time_s_(end) - Ch2_frame_time_sec(end)) < 1
    warning('Analog sync pulse train does not appear to have been stopped before end of recording.')
end

DIO_diff = diff(DATA.DIO0);
DIO_up_ix = find(DIO_diff > 0);
DIO_frame_time_sec = DATA.Time_s_(DIO_up_ix);
DIO_n_frames = length(DIO_frame_time_sec);
DIO_frame_interval_sec = median(diff(DIO_frame_time_sec));
DIO_video_sFreq = 1/DIO_frame_interval_sec;

fprintf('DIO0: %d frames, %1.2f Hz\n',DIO_n_frames, DIO_video_sFreq)

if (DATA.Time_s_(end) - DIO_frame_time_sec(end)) < 1
    warning('DIO sync pulse train does not appear to have been stopped before end of recording.')
end

%% optimized version
v1  = VideoReader(video_file{1});
v2  = VideoReader(video_file{2});
vw = VideoWriter(out_video_file);
open(vw)

% ----- time vector setup
start_time = max(Ch2_frame_time_sec(1), DIO_frame_time_sec(1)) + 1; % only start when both cameras running
end_time = min(Ch2_frame_time_sec(end), DIO_frame_time_sec(end)) - 1;
Ch2_keep = find(Ch2_frame_time_sec >= start_time & Ch2_frame_time_sec <= end_time); % frame idxs to keep
DIO_keep = find(DIO_frame_time_sec >= start_time & DIO_frame_time_sec <= end_time);
if length(Ch2_keep) ~= length(DIO_keep)
    warning('Different number of frames between Ch2 and DIO even after aligning')
end
n_frames = length(Ch2_keep);
frame_time_sec = Ch2_frame_time_sec(Ch2_keep); % combined time vector

% restrict to trials if requested
if cfg.trials_only

    ReadObstacleExpKeys; % reads frame times to keep of Ch2 video
    
    Ch2_trial_keep = zeros(size(Ch2_keep)); DIO_trial_keep = zeros(size(DIO_keep));
    Ch2_trial_number = zeros(size(Ch2_keep));
    for iTrial = 1:length(start_frames)

        Ch2_trial_keep(Ch2_keep >= (start_frames(iTrial) - cfg.trial_padding) & ...
            Ch2_keep <= (end_frames(iTrial) + cfg.trial_padding)) = 1;

        Ch2_trial_number(Ch2_keep >= (start_frames(iTrial) - cfg.trial_padding) & ...
            Ch2_keep <= (end_frames(iTrial) + cfg.trial_padding)) = iTrial;

        DIO_trial_keep(DIO_keep >= (start_frames(iTrial) - cfg.trial_padding) & ...
            DIO_keep <= (end_frames(iTrial) + cfg.trial_padding)) = 1;

    end

end

% advance video readers if needed
Ch2_frame_count = Ch2_keep(1);
for iFrame = 1:Ch2_frame_count-1
    readFrame(v1);
end

DIO_frame_count = DIO_keep(1);
for iFrame = 1:DIO_frame_count-1
    readFrame(v2);
end

% ----- ONE-TIME graphics setup ------------------------------------------
fig        = figure(101); clf(fig)
axVid1     = subplot(4,4,[1 2 5 6], 'Parent', fig);
axVid2     = subplot(4,4,[3 4 7 8], 'Parent', fig);
axPoint    = subplot(4,4, 9:12, 'Parent', fig);
axTrace    = subplot(4,4, 13:16, 'Parent', fig);

% Static trace once
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(axTrace,DATA.Time_s_(GIX),DATA.smoothed(GIX));  axis(axTrace,'tight'); hold(axTrace,'on')
set(axTrace,'XTickLabel',{});

% Pre-create (and remember) the graphics objects we’ll update
hImg1  = image(nan(v1.Height,v1.Width),'Parent',axVid1);     % empty image
hImg2  = image(nan(v2.Height,v2.Width),'Parent',axVid2);      % empty image
hPt2  = plot(axPoint,nan,nan,'b.','MarkerSize',10); hold(axPoint, 'on');
hLinePt = plot(axPoint, [nan nan], [-2 6], 'r-', 'LineWidth', 1); hold(axPoint, 'off');
hPt3 = plot(axTrace, nan, nan, 'r.', 'MarkerSize', 15); hold(axTrace, 'off'); 
title(axPoint,digilent_fname)

axis(axVid1, 'off'); axis(axVid2, 'off');
set(axVid1, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
set(axVid2, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'Box', 'off');
set(axPoint, 'YLim', [-2 6], 'YLimMode', 'manual');
set(axTrace, 'YLim', [-2 6], 'YLimMode', 'manual');


for iFrame = 1:n_frames
    frame_sec = frame_time_sec(iFrame);

    if ~Ch2_trial_number(iFrame)
        readFrame(v1);
        readFrame(v2);
        Ch2_frame_count = Ch2_frame_count + 1;
        DIO_frame_count = DIO_frame_count + 1;
        continue;
    end
    
    vid1Frame  = readFrame(v1);
    hImg1.CData = vid1Frame;
    title(axVid1,sprintf('Frame %d (t = %.4f s, trial %d)',Ch2_frame_count,frame_sec,Ch2_trial_number(iFrame)))

    vid2Frame  = readFrame(v2);
    hImg2.CData = vid2Frame;
    title(axVid2,sprintf('Frame %d (t = %.4f s)',DIO_frame_count,frame_sec))

    rec_ix = find(DATA.Time_s_ >= frame_sec, 1, 'first');

    f1 = max(1, rec_ix - 1000); % first frame of moving window
    f2 = max(1, rec_ix + 1000);
    set(hPt2,'XData',DATA.Time_s_(f1:f2),'YData',DATA.smoothed(f1:f2));
    set(hLinePt, 'XData', [frame_sec frame_sec]);
    set(axPoint, 'XLim', [DATA.Time_s_(f1) DATA.Time_s_(f2)]);

    set(hPt3, 'XData', frame_sec, 'YData', DATA.smoothed(rec_ix));

    drawnow limitrate               % flush graphics efficiently
    frame = getframe(fig);
    writeVideo(vw,frame);           % encode

    Ch2_frame_count = Ch2_frame_count + 1;
    DIO_frame_count = DIO_frame_count + 1;
end

close(vw)

%%
% Plot the data along with the video, frame by frame.

% If there is a DLC file, then load this and make a more fancy plot of the
% data with the movement.

dlc = readtable(dlc_file,'NumHeaderLines',2);
