function [th,back_med] = makeStaticBackgroundModel(varargin)
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input:



%Output:

default_type  = 'files';

default_rect  = [];
default_MSK  = [];
default_SET_TH  = true;
default_limits  = []; % limits for setting TH
default_num_th = 1; % number of th limits (for example - food, body, fish eyes);
default_video = []; % video file
default_dark_on_light = 1; % video file

% parse
vars = inputParser;
addParameter(vars,'Type',default_type);
addParameter(vars,'rect',default_rect);
addParameter(vars,'Msk',default_MSK);
addParameter(vars,'SET_TH',default_SET_TH);
addParameter(vars,'limits',default_limits);
addParameter(vars,'num_th',default_num_th);
addParameter(vars,'vid',default_video);
addParameter(vars,'dark_on_light',default_dark_on_light);



parse(vars,varargin{:})

Type = vars.Results.Type;
rect =  vars.Results.rect;
Msk = vars.Results.Msk;
SET_TH = vars.Results.SET_TH;
limits  =  vars.Results.limits;
num_th = vars.Results.num_th;
vid = vars.Results.vid;
dark_on_light = vars.Results.dark_on_light;

if isempty(limits)
    limits = 1:0.05:2;
end
%...........Local Variable definitions..........


%.................Main Function.................
% try and extract objects from movie
%make a background model:

if strcmp(Type,'files')
    
    files = dir('file*.mat');
    Numf = length(files);
    
    load(files(1).name);
    

    
    [y,x,~] = size(data);
    
    if isempty(rect)
        rect = [0 0 x y];
    end
    
    pic = data(:,:,25);
    pic = imcrop (pic, rect);
    [y,x] = size(pic);
    
    if isempty(Msk)
        Msk = true(size(pic));
    end
    
    Step = round(Numf/100);
    back = uint8(zeros(y,x,length((1:Step:Numf))));
    back(:,:,1) = data(:,:,25);
    j = 2;
    for i = 1:Step:Numf
        name = sprintf('file%d.mat',i);
        load(name);
        data = squeeze(data);
        pic = data(:,:,25);
        pic = imcrop (pic, rect);
        %         pic = rgb2gray(pic);
        pic(~Msk) = 0;
        back(:,:,j) = pic;
        j = j+1;
    end

back_med = mode(single(back),3);
figure(); 
imagesc(back_med); colormap gray; axis image; 
% clear back

elseif strcmp(Type,'movie')
   
    
    pic = rgb2gray(readFrame(vid));
    Numf = vid.NumFrames;
    
    [y,x] = size(pic);
    
    if isempty(rect)
        rect = [0 0 x y];
    end
    
    pic = imcrop (pic, rect);
    [y,x] = size(pic);
    
    if isempty(Msk)
        Msk = true(size(pic));
    end
    
    Step = round(Numf/300);
    back = uint8(zeros(y,x,length((1:Step:Numf))));
%     back(:,:,1) = data(:,:,25);
    j = 1;
    for i = 1:Step:Numf
        vid.CurrentTime = (i/vid.FrameRate);
        pic = rgb2gray(readFrame(vid));
        pic = imcrop (pic, rect);
        %         pic = rgb2gray(pic);
        pic(~Msk) = 0;
        back(:,:,j) = pic;
        j = j+1;
    end

back_med = mode(single(back),3);
figure(); 
imagesc(back_med); colormap gray; axis image; 
end
    
 
%clear very dark spots on the edges:
% back_med(back_med<50) = back_med(back_med<50)/2;

for i = 1: num_th
    SET_TH = 1;
    close all;
    while SET_TH
        ax = [];
        th(i) = setTH(pic,back_med,limits,'dark_on_light',dark_on_light);
        %look at results
        ii = randi(size(back,3),1,5);
        frames = back(:,:,ii);
        for f = 1:size(frames,3)
            pic = frames(:,:,f);
            if dark_on_light
                Ratio = back_med./single(pic);
                Ratio_rev = single(pic)./ back_med;

            else
                Ratio = single(pic)./back_med;
                Ratio_rev = back_med./single(pic);

            end
            
            bin_img = Ratio > th(i);
            
            figure('Units','Normalized','Position',[0.05 0.1 0.8 0.8]);

            ax(1) = subplot(4,2,1);
            imagesc(pic); axis image; colormap gray; axis off;
            title('original image');
            
            ax(2) = subplot(4,2,2);
            imagesc(back_med); axis image; colormap gray; axis off;
            title('background');

            ax(3) = subplot(4,2,3);
            imagesc(Ratio); axis image; colormap gray; axis off;
            title('Ratio');
            
            ax(4) = subplot(4,2,4);
            imagesc(bin_img); axis image; colormap gray; axis off;
            title('binned image');
    
            ax(5) = subplot(4,2,5);
            imagesc(Ratio_rev); axis image; colormap gray; axis off;
            title('reverse ratio');
            ax(6) = subplot(4,2,6);
            imagesc(Ratio_rev> th(i)); axis image; colormap gray; axis off;
            title('reverse binned image');
            
            ax(7) = subplot(4,2,7);
            combined_binned = Ratio_rev> th(i) | bin_img;
            imagesc(combined_binned); axis image; colormap gray; axis off;
            title('combined binned image');
            
            LB = 100; UB = 10000;
            opened_binned = xor(bwareaopen(combined_binned,LB), bwareaopen(combined_binned,UB));
            ax(8) = subplot(4,2,8);
            imagesc(opened_binned); axis image; colormap gray; axis off;
            title('final');
            linkaxes(ax,'xy');

        end
        SET_TH = input('try to set th again? 1 - yes, 0 no - save cuurent th\n');
    end
end

%............Call for local functions...........




