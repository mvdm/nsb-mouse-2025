% Synchronize the video with the mice.
% Cowen 2024
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% define the files.
close all
digilent_file = 'C:\data\nsb2025\M287-2025-07-03_obstacle\acq0020.csv';
video_file = 'C:\data\nsb2025\M287-2025-07-03_obstacle\pilot0001_qtr.mp4';
dlc_file = []; % The deep lab cut csv file if it exists
out_video_dir = "C:\Temp\";
vid_limits_x = [0 960];
vid_limits_y = [0 600];
% dlc_file = 'C:\MBL DATA\test1\recording_testdownsampledDLC_resnet50_hand-testJun30shuffle1_4000.csv';
USE_DIO = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,digilent_fname] = fileparts(digilent_file);
out_video_file = fullfile(out_video_dir,digilent_fname);
% load the digilent data and get the timestamps for each video frame. 
DATA = readtable(digilent_file);
%DATA.smoothed = movmean(DATA.Channel1_V_,2000);
DATA.smoothed = DATA.Channel1_V_;

figure
plot(DATA.Time_s_, DATA.smoothed)


if USE_DIO
    figure;
    plot(DATA.Time_s_,DATA.Channel2_V_)
    d = diff(DATA.Channel2_V_);
else
    figure;
    plot(DATA.Time_s_,DATA.Channel2_V_)
    d = diff(DATA.Channel2_V_ > 1.5);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% Validate and compute sample rates and num frames.
up_ix = find(d>0);
frame_time_sec = DATA.Time_s_(up_ix);
n_frames = length(frame_time_sec);
frame_interval_sec = median(diff(frame_time_sec));
video_sFreq = 1/frame_interval_sec;

fprintf('%d frames, %1.2f Hz\n',n_frames, video_sFreq)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
% load the video file. Check to ensure the number of frames in the video
% data corresponds to the number of frames in the digilent file.
v = VideoReader(video_file);
vw = VideoWriter(out_video_file);
open(vw)

fh = figure(101);
clf
subplot(2,2,3:4)
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(DATA.Time_s_(GIX),DATA.smoothed(GIX))
axis tight
hold on

frame_cnt = 1000;
while hasFrame(v)
    frame_sec = frame_time_sec(frame_cnt);

    subplot(2,2,1)
    % currAxes = axes;

    vidFrame = readFrame(v);
    % image(vidFrame,"Parent",currAxes)
    image(vidFrame)
    xlim(vid_limits_x)
    ylim(vid_limits_y)
    % currAxes.Visible = "off";

    title(sprintf('Frame %d Sec %1.4f',frame_cnt,frame_sec))

    % pause(1/v.FrameRate)

    subplot(2,2,2)
    % find the record closest to this frame.
    rec_ix = find(DATA.Time_s_ >= frame_sec,1,'first');
    plot(frame_sec,DATA.smoothed(rec_ix),'b.')
    %hold on
    title(digilent_fname)

    subplot(2,2,3:4)
    plot(frame_sec,DATA.smoothed(rec_ix),'bo','MarkerSize',3)
    
    frame = getframe(gcf);
    writeVideo(vw,frame)

    frame_cnt = frame_cnt + 1;

    if mod(frame_cnt,200) == 0
        fprintf('%d ',frame_cnt)
    end

end
close(vw)
close(v)

%% optimized version
v  = VideoReader(video_file);
vw = VideoWriter(out_video_file);
open(vw)

% ----- ONE-TIME graphics setup ------------------------------------------
fig        = figure(101); clf(fig)
axVid      = subplot(2,2,1,  'Parent',fig);
axPoint    = subplot(2,2,2,  'Parent',fig);
axTrace    = subplot(2,2,3:4,'Parent',fig);

% Static trace once
GIX = DATA.Time_s_ >= frame_time_sec(1) & DATA.Time_s_ <= frame_time_sec(end);
plot(axTrace,DATA.Time_s_(GIX),DATA.smoothed(GIX));  axis(axTrace,'tight'); hold(axTrace,'on')

% Pre-create (and remember) the graphics objects we’ll update
hImg  = image(nan(v.Height,v.Width),'Parent',axVid);      % empty image
hPt2  = plot(axPoint,nan,nan,'b.','MarkerSize',11);  hold(axPoint,'on')
hPt3  = plot(axTrace,nan,nan,'bo','MarkerSize',4);

frame_cnt = 1;
while hasFrame(v)
    frame_sec = frame_time_sec(frame_cnt);
    vidFrame  = readFrame(v);

    % --- update graphics objects, don’t create new ones ---
    hImg.CData = vidFrame;
    title(axVid,sprintf('Frame %d  (t = %.4f s)',frame_cnt,frame_sec))

    rec_ix = find(DATA.Time_s_ >= frame_sec, 1, 'first');

    f1 = max(1, rec_ix-500); % first frame of moving window
    set(hPt2,'XData',DATA.Time_s_(f1:rec_ix),'YData',DATA.smoothed(f1:rec_ix));

    set(hPt3,'XData',[hPt3.XData, frame_sec], ...
             'YData',[hPt3.YData, DATA.smoothed(rec_ix)]);

    title(axPoint,digilent_fname)

    drawnow limitrate               % flush graphics efficiently
    frame = getframe(fig);
    writeVideo(vw,frame);           % encode

    frame_cnt = frame_cnt + 1;
end

close(vw)

%%
% Plot the data along with the video, frame by frame.

% If there is a DLC file, then load this and make a more fancy plot of the
% data with the movement.

dlc = readtable(dlc_file,'NumHeaderLines',2);
