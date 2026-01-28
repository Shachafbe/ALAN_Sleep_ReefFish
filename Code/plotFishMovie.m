function  plotFishMovie(coor,mov,varargin)
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input: 



%Output: 

default_fish_i  = 1;
default_points_i  = 1:2;

default_start_f = 1;
default_end_f  = 1:10000;
default_data_per_point  = 3; % amount of info per point (x,y,quality)
default_points_per_fish  = 6; % number of points per fish



% parse
vars = inputParser;
addParameter(vars,'fish_i',default_fish_i);
addParameter(vars,'points_i',default_points_i);
addParameter(vars,'start_f',default_start_f);
addParameter(vars,'end_f',default_end_f);
addParameter(vars,'data_per_point',default_data_per_point);
addParameter(vars,'points_per_fish',default_points_per_fish);


parse(vars,varargin{:})

fish_i = vars.Results.fish_i;
points_i = vars.Results.points_i;
start_f = vars.Results.start_f;
end_f = vars.Results.end_f;
data_per_point = vars.Results.data_per_point;
points_per_fish = vars.Results.points_per_fish;


%...........Local Variable definitions..........
Cmap = LineMap;

%.................Main Function.................


% create movie objects


% if MAKEMOVIE
%         movie_name = ['compare_light_dark_and_speed_fish1_smoothed'];
%         obj = VideoWriter([movie_name],'MPEG-4');
%         obj.FrameRate = 25;
%         obj.Quality = 100;
%         open(obj);
% end

H = figure('Units','Normalized','Position',[0.2 0.2 0.6 0.6]);
pause(0.1);
% Tail = 50;
% Ylim = [0 25];
mov.CurrentTime = start_f/mov.FrameRate;

for i = start_f:end_f
    G = readFrame(mov);
%     subplot(2,2,1);
    if ishghandle(H)

        imagesc(G); colormap gray; axis image; hold on;
%     else
%        disp('Break by user');
%        return
%     end

    for f = 1:length(fish_i) % loop over fish
        fi = fish_i(f);
        
        % take all points of the current fish:
        curr_points = squeeze(coor(i,f,points_i,:));
        x = curr_points(:,1);
        y = curr_points(:,2);
        

        plot(x,y,'.','Color',Cmap(fi,:));

       
    end
    title(['frame: ',num2str(i)]);
    pause(0.01);
    hold off
    else 
        retuen
    end

end

        
            
%             
%         
%    
% 
% if MAKEMOVIE
%    close(obj); 
% end




%............Call for local functions...........




