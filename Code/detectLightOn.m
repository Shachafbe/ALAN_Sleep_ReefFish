function [on_off, pix_values] = detectLightOn(filename)
% gets the name of an image file, interactively draws a region of interst
% and detets all frames when light is on


% define th
TH = 1.05; % anything that is 5% obove the baseling will be set as light on


% create a video object
obj = VideoReader(filename);   

% get firat frame
frame = readFrame(obj);
frame = rgb2gray(frame);

% create a region of interset by drawing a ploygon:
ROI = roipoly(frame);

estimated_nframes = ceil(obj.Duration*obj.FrameRate);
pix_values = zeros(1,estimated_nframes);

% add first avg pixel value in region
pix_values(1) = mean(frame(ROI));
% all_frames = zeros(size(frame,1),size(frame,2),estimated_nframes,'uint8');
% all_frames(:,:,1) = frame;
% loop over all frames and  extract pixel values in that region:
p = 2;
while hasFrame(obj)
    frame = readFrame(obj);
    frame = rgb2gray(frame);
%     all_frames(:,:,p) = frame;
    pix_values(p) = mean(frame(ROI));
    p = p+1;
end

if p<length(pix_values)
    pix_values(p+1:end) = [];
end


% define baseline:
baseline = mode(pix_values(1:100:end));

% fine when avg pixels are brighter then the mode
on_off  = pix_values > TH*baseline;
figure();
plot(pix_values);
hold on
plot(find(on_off),pix_values(on_off),'k.');
xlabel('frames'); ylabel('intenisty');
box off


