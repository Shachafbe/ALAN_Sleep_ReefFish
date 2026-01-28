function [xc,yc,Sc,Ang,q] = preProcessSleepingFisData(coor,fish_id, points_on_fish,horz_axis,varargin)
%description: runs over all coordinates, calculates speed and bidy
%postures. remoevs (interpulates) over points_on_fish that seems like outliers


%Input: coor - cell array of fish coordinates, within each cell a matrix of
%              T x dims (num points_on_fish x 3 which is x,y,quality per point)
%              fish_id -  matrix of size(coor) with integers corresponsing
%              to fish id.
%       points_on_fish - vector of 1xp which are the indices of points_on_fish on the fish
%               p1 = 1; % head point, p2 = 2; % back point, p3 = 3; % caudal fin

%       horz_axis - matrix of numfish x 2 calculated vectors (per fish) representing horz axis

% Optional inputs - speed_th, quality_th, filt_win (see defaults below)



%Output: 



%...........Local Variable definitions..........

% set defaults
default_speed_th =  40; % speed th to remove
default_quality_th = 0.9; % quality th
default_filt_win = 10; % filtering window
defaultPLOT = 0;
defaultSAVE = 0;




% parse
vars = inputParser;
addParameter(vars,'speed_th',default_speed_th);
addParameter(vars,'quality_th',default_quality_th);
addParameter(vars,'filt_win',default_filt_win);
addParameter(vars,'PLOT',defaultPLOT);
addParameter(vars,'SAVE',defaultSAVE);



parse(vars,varargin{:})
% make variables
speed_th =  vars.Results.speed_th; %m
quality_th =  vars.Results.quality_th; %m
filt_win =  vars.Results.filt_win; %m
PLOT = vars.Results.PLOT;
SAVE = vars.Results.SAVE;



Num_points_per_fish = length(points_on_fish);



%.................Main Function.................

%% Try to clean up the data a little bit:

% extract head coordinates and quality ratings from the data matrices:

% light data
x = cell(size(coor)); % head x coor of ligth files
y = cell(size(coor)); % y coor
q = cell(size(coor)); % quality of tracking

S = cell(size(coor)); % variable for calculated speed
Ang = cell(size(coor)); % variable for body angles:



for i = 1:length(x(:))
    % data is [eye, dorsal fish, tail fork]
    if isempty(coor{i})
        continue
    end
        
    x{i} = coor{i}(:,1:3:end); % extract only relevant columns all points_on_fish per fish
    y{i} = coor{i}(:,2:3:end);
    q{i} = coor{i}(:,3:3:end);
    
    S{i} = calculateNorm(cat(3,x{i}(3:end,:)-x{i}(1:end-2,:),...
        y{i}(3:end,:)-y{i}(1:end-2,:)))/2;
    S{i} = [zeros(1,Num_points_per_fish); S{i}; zeros(1,Num_points_per_fish)];
    
    % calculate the angles from horizontal axis: [1 - head - dorsal,

    
    % points_on_fish on the fish: p1 = 1; % head point, p2 = 2; % back point, p3 = 3; % caudal fin
    
    
   % var for angles:
    angs_of_body = zeros(size(x{i},1),Num_points_per_fish);
        
    % define main body vevtor: p3-p1 (head to tail)
    vec_body = [x{i}(:,points_on_fish(3)) y{i}(:,points_on_fish(3))] - [x{i}(:,points_on_fish(1)) y{i}(:,points_on_fish(1))];

    % find instances where fish is pointing right and not left
    jj = vec_body(:,1) < 0;
    
    % define additional vectors:
    vec_back_tail = [x{i}(:,points_on_fish(3)) y{i}(:,points_on_fish(3))] - [x{i}(:,points_on_fish(2)) y{i}(:,points_on_fish(2))];
    vec_eye_back = [x{i}(:,points_on_fish(2)) y{i}(:,points_on_fish(2))] - [x{i}(:,points_on_fish(1)) y{i}(:,points_on_fish(1))];
    
    % for angle calculations switch direction so the fish is always
    % pointing left (otherwise we wont get consitent angles)
    vec_body(jj,1) = -(vec_body(jj,1));
    vec_back_tail(jj,1) = -(vec_back_tail(jj,1));
    vec_eye_back(jj,1) = -(vec_eye_back(jj,1));
    %     vec_top_bottom_c(:,1) = abs(vec_top_bottom_c(:,1));
    
    % check if this is a fish on the left or on the right for using the
    % correct reference:
    ref_x = horz_axis(fish_id(i),:);

    
    % calculate angles:
    % body angle from horizontal
    [~,angs_of_body(:,1)] = angOfVectors(ref_x, vec_body,1);    
    % back_caudal top from body
    [~,angs_of_body(:,2)] = angOfVectors(vec_body, vec_back_tail, 1);
    % back_caudal bottom from body
    [~,angs_of_body(:,3)] = angOfVectors(vec_body, vec_eye_back, 1);
    % caudal bottom-top from body
    [~,angs_of_body(:,4)] = angOfVectors(vec_eye_back, vec_back_tail, 1);
    
    % make all angs (clockwise from ref)between -180 and 180 
    angs_of_body(angs_of_body>180) = angs_of_body(angs_of_body > 180) - 360;
    
    % save angles to var
    Ang{i} = angs_of_body;

end

%%

% decide speed and quality thresholds to interpulate over:


disp([' removing points_on_fish with speed > ',num2str(speed_th),...
    ' and quality < ',num2str(quality_th)]);

% clean up the data by interpulating over spurious points_on_fish:

% variables for interpulated data
xc = cell(size(coor)); % x for cleaned up light data
yc = cell(size(coor)); % y
Sc = cell(size(coor)); % speed


% loop over experiments and then over fish:

% for light experiments:
for j = 1:length(x(:)) % all experiments
    if isempty(x{j})
        continue
    end
    for i = 1:Num_points_per_fish

        x_temp = x{j}(:,i); y_temp = y{j}(:,i);
        
        % find suspicious points_on_fish to correct
        jj =  find(S{j}(:,i) > speed_th | q{j}(:,i) < quality_th | x_temp<0);

        
        % start interpulating before and after the points_on_fish:
        temp_qs = false(size(x));
        for pp = -2:2
            ii = jj+pp;
            ii(ii<1 | ii > length(x)) = [];
            temp_qs(ii) = true;
        end
        
        % interpulate over suspicious points_on_fish:
        x_temp(temp_qs) = interp1(find(~temp_qs),x_temp(~temp_qs),find(temp_qs));
        y_temp(temp_qs) = interp1(find(~temp_qs),y_temp(~temp_qs),find(temp_qs));
        
        % smooth the data 
        xc{j}(:,i) = smoothdata(x_temp,'sgolay',filt_win);
        yc{j}(:,i) = smoothdata(y_temp,'sgolay',filt_win);
        
%         xcl{j}(:,i) = x;
%         ycl{j}(:,i) = y;
    end
    % calculate speed estimation of the corrected data:
    Sc{j} = calculateNorm(cat(3,xc{j}(3:end,:)-xc{j}(1:end-2,:),...
        yc{j}(3:end,:)-yc{j}(1:end-2,:)));
    Sc{j} = [zeros(1,Num_points_per_fish);Sc{j};zeros(1,Num_points_per_fish)];
    disp(['Finished correcting experiment #',num2str(j)]);
end


%............Call for local functions...........




