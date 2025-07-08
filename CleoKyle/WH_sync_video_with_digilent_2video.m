% Synchronize the video with the mice.
% Cowen 2024
% MvdM 2025 adapted for 2 video streams
% analog ch2: video stream 1
% DIO 0: video stream 2

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% define the files.
close all
digilent_file = 'C:\data\nsb2025\M284-2025-07-07_obstacle\acq0020.csv';
video_file{1} = 'C:\data\nsb2025\M284-2025-07-07_obstacle\video1_qtr.mp4';
video_file{2} = 'C:\data\nsb2025\M284-2025-07-07_obstacle\video2_qtr.mp4';
dlc_file = []; % The deep lab cut csv file if it exists
out_video_dir = "C:\Temp\";
vid_limits_x = [0 960];
vid_limits_y = [0 600];
% dlc_file = 'C:\MBL DATA\test1\recording_testdownsampledDLC_resnet50_hand-testJun30shuffle1_4000.csv';
% USE_DIO = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,digilent_fname] = fileparts(digilent_file);
out_video_file = fullfile(out_video_dir,digilent_fname);

% load the digilent data and get the timestamps for each video frame. 
DATA = readtable(digilent_file);
DATA.smoothed = movmean(DATA.Channel1_V_,50); % this is the fiber signal
%DATA.smoothed = DATA.Channel1_V_;

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
start_time = max(Ch2_frame_time_sec(1), DIO_frame_time_sec(1)); % only start when both cameras running
end_time = min(Ch2_frame_time_sec(end), DIO_frame_time_sec(end));
Ch2_keep = find(Ch2_frame_time_sec >= start_time & Ch2_frame_time_sec <= end_time); % frame idxs to keep
DIO_keep = find(DIO_frame_time_sec >= start_time & DIO_frame_time_sec <= end_time);
if length(Ch2_keep) ~= length(DIO_keep)
    error('Different number of frames between Ch2 and DIO even after aligning')
end
n_frames = length(Ch2_keep);
frame_time_sec = Ch2_frame_time_sec(Ch2_keep); % combined time vector

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
axVid1     = subplot(2,2,1, 'Parent', fig);
axVid2     = subplot(2,2,2, 'Parent', fig);
axPoint    = subplot(2,2,3, 'Parent', fig);
axTrace    = subplot(2,2,4,'Parent', fig);

% Static trace once
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(axTrace,DATA.Time_s_(GIX),DATA.smoothed(GIX));  axis(axTrace,'tight'); hold(axTrace,'on')

% Pre-create (and remember) the graphics objects we’ll update
hImg1  = image(nan(v1.Height,v1.Width),'Parent',axVid1);      % empty image
hImg2  = image(nan(v2.Height,v2.Width),'Parent',axVid2);      % empty image
hPt2  = plot(axPoint,nan,nan,'b.','MarkerSize',11);  hold(axPoint,'on')
hPt3  = plot(axTrace,nan,nan,'bo','MarkerSize',4);

for iFrame = 1:n_frames
    frame_sec = frame_time_sec(iFrame);
    
    vid1Frame  = readFrame(v1);
    hImg1.CData = vid1Frame;
    title(axVid1,sprintf('Frame %d (t = %.4f s)',Ch2_frame_count,frame_sec))

    vid2Frame  = readFrame(v2);
    hImg2.CData = vid2Frame;
    title(axVid2,sprintf('Frame %d (t = %.4f s)',DIO_frame_count,frame_sec))

    rec_ix = find(DATA.Time_s_ >= frame_sec, 1, 'first');

    f1 = max(1, rec_ix-500); % first frame of moving window
    set(hPt2,'XData',DATA.Time_s_(f1:rec_ix),'YData',DATA.smoothed(f1:rec_ix));

    set(hPt3,'XData',[hPt3.XData, frame_sec], ...
             'YData',[hPt3.YData, DATA.smoothed(rec_ix)]);

    title(axPoint,digilent_fname)

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
