function new_coor = rearrangeCoorData(coor,varargin)
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input: 



%Output: 

default_num_fish  = 6; 
default_data_per_point  = 3; % amount of info per point (x,y,quality)
default_points_per_fish  = 6; % number of points per fish



%parse
vars = inputParser;
addParameter(vars,'num_fish',default_num_fish);
addParameter(vars,'data_per_point',default_data_per_point);
addParameter(vars,'points_per_fish',default_points_per_fish);


parse(vars,varargin{:})

num_fish = vars.Results.num_fish;
data_per_point = vars.Results.data_per_point;
points_per_fish = vars.Results.points_per_fish;


%...........Local Variable definitions..........
[T,cols] = size(coor);

new_coor = zeros(T,num_fish,points_per_fish,data_per_point);

pp = 1;
for n = 1:num_fish
    for p = 1:points_per_fish
        for d = 1:data_per_point
            new_coor(:,n,p,d) = coor(:,pp);
            pp = pp+1;
        end
    end
end



%.................Main Function.................



%............Call for local functions...........




