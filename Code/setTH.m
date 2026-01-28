function th = setTH(pic,back_med,limits,varargin)



default_dark_on_light = 1; % dark fish on light back ground or the other way arounf



% parse
vars = inputParser;
addParameter(vars,'dark_on_light',default_dark_on_light);
parse(vars,varargin{:})

dark_on_light = vars.Results.dark_on_light;

%loop to find empirical TH
if nargin < 3
   limits = (1:0.05:3);
end
    
thresh = limits;
ObjN = zeros(1,size(thresh,2));

figure('Units','Normalized','Position',[0.2 0.1 0.6 0.8])

for t = 1:size(thresh,2)
    if dark_on_light
        bin_img = (back_med./single(pic)>thresh(t));
        reverse_img = (single(pic)./back_med>thresh(t));
    else
         bin_img = (single(pic)./back_med>thresh(t));
         reverse_img = (back_med./single(pic))>thresh(t);

    end
        bin_img = bwareaopen(bin_img,5);
        reverse_img = bwareaopen(reverse_img,5);

        %calculate statistics of objects
%         STATS = regionprops(bin_img,'PixelIdxList','Area','Centroid','MajorAxisLength'...
%             ,'MinorAxisLength','Orientation','Image','PixelList','BoundingBox');
        
        subplot(1,2,1);
        imagesc(bin_img);colormap gray;
        axis image;
        title(num2str(thresh(t)));
        subplot(1,2,2);
        imagesc(reverse_img);colormap gray;
        axis image;
        title('reverse light polarity');
        pause(); hold off;

%         ObjN(t) = size(STATS,1);
end

% figure();
% plot(thresh,ObjN);
th = input('please enter th: ');
close all

