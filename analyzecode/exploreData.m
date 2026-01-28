% routine to run perliminary analysis of the fish movies 

% define source folder where data files (.csv) and movie files (.avi) are
% at. IMPORTANT: movie files and data files should have the same exact
% names, with only different extentions.

% define colormaps for lines 

Cmap = LineMap;

% folder where .csv and .mp4 files are places
% change to your folder
Folder = 'D:\Sleeping fish\4fish\'; 

% put here the folder where you want files to be saved to (right now only
% .fig files are saved)
save_folder = 'C:\Users\Roy h\Dropbox\Collaborations\CoralFish\For Shachaf\Results\';

% important decide if to save files and make a movie as default (I usually
% leave it as 0, and change it only if I need to).
SAVE = 0;
MAKEMOVIE = 0; % decide if to make a movie or not

cd(Folder); % change to folder location


% put in file names (without extensions):
filenames_light = {'LN1DeepCut_resnet50_lightdarkNov13shuffle1_650000',...
                'LN2DeepCut_~1930-0500',...
                'LN3DeepCut_resnet50_lightdarkNov13shuffle1_650000'};
filenames_dark = {'DN1DeepCut_~1930-0540',...
    'DN2DeepCut_resnet50_lightdarkNov13shuffle1_350000'};


% put in file names (without extensions):
filenames_light = dir('LN*.csv'); % all files starting with LN and a .csv extension
filenames_light = {filenames_light.name};
filenames_dark = dir('DN*.csv'); % all files starting with DN and a .csv extension
filenames_dark = {filenames_dark.name};

orderl = [1 2 3]; % days of the experiment
orderd = [4 5]; 

numl = length(filenames_light); % number of light experiments
numd = length(filenames_dark);  % number of dark experiments

% lood data files and create movie objects:
coorl = cell(1,numl); 
coord = cell(1,numd);

% create also links to video files
vid_l = cell(1,numl); 
vid_d = cell(1,numd);

% loop to load data files and video objects (takes some time)
for i = 1:numl
    coorl{i} = csvread([Folder,filenames_light{i}],3,1);
    if exist([Folder,filenames_light{i}(1:end-4),'.mp4'],'file')
    vid_l{i} = VideoReader([Folder,filenames_light{i}(1:end-4),'.mp4']);
    end
    disp(['Finished light file #',num2str(i)]);
end


for i = 1:numd
    coord{i} = csvread([Folder,filenames_dark{i}],3,1);
    if exist([Folder,filenames_dark{i}(1:end-4),'.mp4'],'file')
        vid_d{i} = VideoReader([Folder,filenames_dark{i}(1:end-4),'.mp4']);
    end
    disp(['Finished dark file #',num2str(i)]);
end
    
% titles for plotting
Titles_light = {'LN1','LN2','LN3'};
Titles_dark = {'DN1','DN2'};


Nf = 4; % number of fish (no including the coordinated used for measuring in the files)
Num_dots = 1; % number of markers on a fish
% Tl = size(coorl,1);
% Td = size(coord,1);
output_per_dot = 3; % number of outputs from DLC per dot (x.y.q) just in case this ever changes

%% Try to clean up the data a little bit:

% extract head coordinates and quality ratings from the data matrices:

% light data
xhl = cell(1,numl); % head x coor of ligth files
yhl = cell(1,numl); % y coor
qhl = cell(1,numl); % quality of tracking

Shl = cell(1,numl); % variable for calculated speed

for i = 1:numl
    xhl{i} = coorl{i}(:,1:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf); % extract only relevant columns
    yhl{i} = coorl{i}(:,2:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf+1);
    qhl{i} = coorl{i}(:,3:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf+2);
    
    Shl{i} = calculateNorm(cat(3,xhl{i}(3:end,:)-xhl{i}(1:end-2,:),...
        yhl{i}(3:end,:)-yhl{i}(1:end-2,:)));
    Shl{i} = [zeros(1,Nf); Shl{i}; zeros(1,Nf)];
    
end


% dark data
xhd = cell(1,numd); % head x coor of ligth files
yhd = cell(1,numd); % y coor
qhd = cell(1,numd); % quality of tracking

Shd = cell(1,numl); % variable for calculated speed

for i = 1:numd
    xhd{i} = coord{i}(:,1:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf); % extract only relevant columns
    yhd{i} = coord{i}(:,2:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf+1);
    qhd{i} = coord{i}(:,3:output_per_dot*Num_dots:output_per_dot*Num_dots*Nf+2);
    
    Shd{i} = calculateNorm(cat(3,xhd{i}(3:end,:)-xhd{i}(1:end-2,:),...
        yhd{i}(3:end,:)-yhd{i}(1:end-2,:)));
    Shd{i} = [zeros(1,Nf); Shd{i}; zeros(1,Nf)];
end


% decide speed and waulity thresholds to interpulate over:
%(thesse can be manually changed or later decided upon in a more
%informative way/smarter way)
speed_th = 50; % above this speed a point is considered spurious
quality_th = 0.8; % below this quality a point is considered spurious
filt_win = 10; % number of frames to smooth over 

disp([' removing points with speed > ',num2str(speed_th),...
    ' and quality < ',num2str(quality_th)]);

% clean up the data by interpulating over spurious points:

% variables foe interpulated data
xcl = cell(1,numl); % x for cleaned up light data
ycl = cell(1,numl); % y
Scl = cell(1,numl); % speed

xcd = cell(1,numd); % x for clearned up dark data
ycd = cell(1,numd); % y
Scd = cell(1,numd); % speed

% loop over experiments and then over fish:

% for light experiments:
for j = 1:numl % all experiments
    for i = 1:Nf
        
        x = xhl{j}(:,i); y = yhl{j}(:,i);
        
        % find suspicious points to correct
        jj =  find(Shl{j}(:,i) > speed_th | qhl{j}(:,i) < quality_th);
        
        % start interpulating before and after the points:
        temp_qs = false(size(x));
        for pp = -2:2
            ii = jj+pp;
            ii(ii<1 | ii > length(x)) = [];
            temp_qs(ii) = true;
        end
        
        % interpulate over suspicious points:
        x(temp_qs) = interp1(find(~temp_qs),x(~temp_qs),find(temp_qs));
        y(temp_qs) = interp1(find(~temp_qs),y(~temp_qs),find(temp_qs));
        
        % smooth the data 
        xcl{j}(:,i) = smoothdata(x,'sgolay',filt_win);
        ycl{j}(:,i) = smoothdata(y,'sgolay',filt_win);
        
    end
    % calculate speed estimation of the corrected data:
    Scl{j} = calculateNorm(cat(3,xcl{j}(3:end,:)-xcl{j}(1:end-2,:),...
        ycl{j}(3:end,:)-ycl{j}(1:end-2,:)));
    Scl{j} = [zeros(1,Nf);Scl{j};zeros(1,Nf)];
    disp(['Finished correcting light experiment #',num2str(j)]);
end

% for dark experiments:
for j = 1:numd % all experiments
    for i = 1:Nf
        
        x = xhd{j}(:,i); y = yhd{j}(:,i);
        
        % find suspicious points to correct
        jj =  find(Shd{j}(:,i) > speed_th | qhd{j}(:,i) < quality_th);
        
        % start interpulating before and after the points:
        temp_qs = false(size(x));
        for pp = -2:2
            ii = jj+pp;
            ii(ii<1 | ii > length(x)) = [];
            temp_qs(ii) = true;
        end
        
        % interpulate over suspicious points:
        x(temp_qs) = interp1(find(~temp_qs),x(~temp_qs),find(temp_qs));
        y(temp_qs) = interp1(find(~temp_qs),y(~temp_qs),find(temp_qs));
        
        % smooth the data 
        xcd{j}(:,i) = smoothdata(x,'sgolay',filt_win);
        ycd{j}(:,i) = smoothdata(y,'sgolay',filt_win);
        
    end
    
    % calculate speed estimation of the corrected data:
    Scd{j} = calculateNorm(cat(3,xcd{j}(3:end,:)-xcd{j}(1:end-2,:),...
        ycd{j}(3:end,:)-ycd{j}(1:end-2,:)));
    Scd{j} = [zeros(1,Nf);Scd{j};zeros(1,Nf)];
    disp(['Finished correcting dark experiment #',num2str(j)]);
end


%% plot trajectories and speed distributions:

% choose which groups to plot:
groups_l = [1]; % out of 3 groups (i.e. [1,2,3] will plot all 3 light groups, and [] will not plot any)
groups_d = [1]; % out of two groups 

% shades for histogram
shades_l = linspace(1,0.2,numl);
shades_d = linspace(1,0.2,numd);

num_rows = max([length(groups_d),length(groups_l)]); % number of columns in subplots
num_cols = 2; % always 2 since we compare dark and light

th_quality_for_traj_plot = 0.9;
% position map:
% start ploting light night groups (from the left)
f1 = figure();
p = 1;
for l = groups_l % all groups
    subplot(num_rows,num_cols,(p-1)*2+1);
    if ~isempty(vid_l{l})
    img = readFrame(vid_l{l});
    imshow(img); 
    end
    hold on;
    for f = 1:Nf
        % remove any spurious data that might clutter  the plots (just for visualizations):
%         
        qsl = Scl{l}(:,f) < 12 & qhl{l}(:,f) > th_quality_for_traj_plot;
        plot(xcl{l}(qsl,f),ycl{l}(qsl,f),'.','MarkerSize',1);

    end
    title(Titles_light{l});
    p = p+1;
end

% add dark night groups (from the right)
p = 1;
for d = groups_d % all groups
    subplot(num_rows,num_cols,(p-1)*2+2);
    if ~isempty(vid_d{d})
    img = readFrame(vid_d{d});
    imshow(img); 
    end
    hold on;
    for f = 1:Nf
        qsd = Scd{d}(:,f) < 12 & qhd{d}(:,f) > th_quality_for_traj_plot;
        plot(xcd{d}(qsd,f),ycd{d}(qsd,f),'.','MarkerSize',1);

    end
    title(Titles_dark{d});
    p = p+1;
end


if SAVE
    name = ['trajectories_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved group trajectories');
end    

% now plot the distributions of speeds:

f2 = figure();

% first light distributions:
p = 1;
for l = groups_l
    for f = 1:Nf
        
        subplot(1,Nf,f);
        histogram(log(Scl{l}(:,f)),-7:0.1:7,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(1,:)*shades_l(p)); hold on;
        box off;
        xlabel('log(Speed)');
        title(['Fish ',num2str(f)]);
    end
    p = p+1;
end

% now dark nigth distributions
p = 1;
for d = groups_d
    for f = 1:Nf
        
        subplot(1,Nf,f);
        histogram(log(Scd{d}(:,f)),-7:0.1:7,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(2,:)*shades_d(p)); hold on;
        box off;
        xlabel('log(Speed)');
    end
    p = p+1;
end

% add legend to last subplot
legend([Titles_light(groups_l), Titles_dark(groups_d)]);
legend boxoff;
%     subplot(2,3,i);
        %     plot(Vorig); hold on;
        %     plot(V);
        %     plot(V); hold on;
        %     histogram(V);

if SAVE
    name = ['speed_distributions_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved speed distributions');
end        


%% look at times of 'up' and 'down' states: 

% define thresholds for "up" and "down" states:
% th = exp([0.01 0.01 0.01]);
th = exp(ones(1,Nf)*0.1);

% for simplicity we use the same th for all fish in all experiments. but we
% can make a smarter desicion according to the speed distributions of the
% experiments we are comparing

% variables for up down states
up_down_l = cell(1,numl);
up_down_d = cell(1,numd);

% variables for start and end of up down of events (time frames), and for
% there length (end-start+1):
lims_l = cell(1,numl);
length_l = cell(1,numl);

lims_d = cell(1,numd);
length_d = cell(1,numd);

% variables for length of low quality events (control analysis): 
length_q_l = cell(3,1); % for length of non quality events
length_q_d = cell(3,1); % for length of non quality events

% define short events to ignore:
short_events = 1; % in frames

% light experiments:
for l = 1:numl
    % temporary varibles per experiment
    temp_up_down = false(size(Scl{l}));
    temp_length = cell(1,3); % for 3 fish
    temp_lims = cell(1,3);
    
    for f = 1:Nf % loop over fish
        temp = Scl{l}(:,f); % take speed
        up_down = temp > th(f); % find up and down
        [st,nd,temp_up_down(:,f)] = findNonZeroSeq(up_down, short_events); % find non zero sequences (remove short events)
        temp_length{f} = nd-st+1; % calcualte their length
        temp_lims{f} = [st;nd];
        
    end
    % save all fish data in experiment variables:
    up_down_l{l} = temp_up_down;
    length_l{l} = temp_length;
    lims_l{l} = temp_lims;
    
end


% dark experiments:
for d = 1:numd
    % temporary varibles per experiment
    temp_up_down = false(size(Scd{d}));
    temp_length = cell(1,3); % for 3 fish
    temp_lims = cell(1,3);
    
    for f = 1:Nf % loop over fish
        temp = Scd{d}(:,f); % take speed
        up_down = temp > th(f); % find up and down
        [st,nd,temp_up_down(:,f)] = findNonZeroSeq(up_down, short_events); % find non zero sequences (remove short events)
        temp_length{f} = nd-st+1; % calcualte their length
        temp_lims{f} = [st;nd];
        
    end
    % save all fish data in experiment variables:
    up_down_d{d} = temp_up_down;
    length_d{d} = temp_length;
    lims_d{d} = temp_lims;
end


%% plot the distributions of up and down states

figure();

% redefine groups to plot:
% choose which groups to plot:
groups_l = [1]; % out of 3 groups (i.e. [1,2,3] will plot all 3 light groups, and [] will not plot any)
groups_d = [1]; % out of two groups 

% shades for histogram
shades_l = linspace(1,0.2,numl);
shades_d = linspace(1,0.2,numd);

p = 1;
for l = groups_l
    for f = 1:Nf
        
        subplot(1,Nf,f);
        histogram(length_l{l}{f},0:1:100,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(1,:)*shades_l(p));
        hold on;
        box off;
        xlabel('Length of activity events [frames]');

    end
end


p = 1;
for d = groups_d
    for f = 1:Nf
        
        subplot(1,Nf,f);
        histogram(length_d{d}{f},0:1:100,'Normalization','Probability',...
            'EdgeColor','none','FaceColor',Cmap(2,:)*shades_d(p));
        hold on;
        box off;
        xlabel('Length of activity events [frames]');
    end
end

% add legend to last subplot
legend([Titles_light(groups_l), Titles_dark(groups_d)]);
legend boxoff;

if SAVE
    name = ['active_in_active_dist_',[Titles_light{groups_l}],...
        [Titles_dark{groups_d}]];
    saveas(gcf,[save_folder,name]);
    disp('saved active and inactive distributions');
end

%% plot fish heatmaps

% for light groups:
simple_maps = LineMap;
nbins = 50; % number of bins in 2d histogram
smooth_val = 0.5; % smoothing factor for 2D histogram

for l = groups_l
    if ~isempty(vid_l{l})
        img = readFrame(vid_l{l});
    else
        img = [];
    end
    figure('Units','Normalized','Position',[0.0 0.0 1 1]);

    for f = 1:Nf
        tempx = xcl{l}(:,f); % get current coordinates:
        tempy = ycl{l}(:,f);
        
        % remove nans
        tempx = tempx(~isnan(tempx));
        tempy = tempy(~isnan(tempy));

        
        subplot(Nf,2,f*2-1);
        
        if ~isempty(img)
            imshow(img);
        end
        hold on; % draw current image
%         scatter(tempx(1:5:end),tempy(1:5:end),0.1,'MarkerEdgeColor','none',...
%             'MarkerFaceColor',Cmap(f,:),'MarkerFaceAlpha',0.1);
        
        % calculate 2D smooth histogram
        [ctr1,ctr2,Counts] = smoothhist2D([tempx tempy],smooth_val,[nbins nbins],0,'none'); 
%        
        [Nrow,Ncol,~] = size(Counts); % number of rows and columns of histogram       
        % create a color mask same size as histogram
        Color = cat(3,repmat(simple_maps(f,1),Nrow,Ncol),repmat(simple_maps(f,2),Nrow,Ncol),...
            repmat(simple_maps(f,3),Nrow,Ncol));
        % draw the color mask
        hh = imagesc(ctr1,ctr2,Color);colormap(CmapRed);
        % change the transperancy of color mask according to 2d hist such
        % that we will see colors only when heat map is high:
        set(hh,'AlphaData',Counts);
        
        title([Titles_light{l},' fish: ',num2str(f)]);
        % draw only a heat map next to the image for better view
        subplot(Nf,2,f*2); 
        hh = imagesc(ctr1,ctr2,Counts);colormap(CmapRed); hold on;
        set(gca,'Xtick',[],'XtickLabel',[],'Ytick',[],'YtickLabel',[]);
        axis image;
        if ~isempty(img)
        axis([0 size(img,2) 0 size(img,1)]); 
        end
        
    end
    if SAVE
       name = ['Heatmaps_',Titles_light{l}];
       saveas(gcf,[save_folder,name]);
       disp(['saved heatmap for ',Titles_light{l}]);
    end
end

% for dark groups:


for d = groups_d
    if ~isempty(vid_d{d})
        img = readFrame(vid_d{d});
    else
        img = [];
    end
    figure('Units','Normalized','Position',[0.0 0.0 1 1]);
    for f = 1:Nf
        tempx = xcd{d}(:,f); % get current coordinates:
        tempy = ycd{d}(:,f);
        
        tempx = tempx(~isnan(tempx)); % get rid of nans
        tempy = tempy(~isnan(tempy));
        
        subplot(Nf,2,f*2-1);
        
        if ~isempty(img)
            imshow(img); 
        end
        hold on; % draw current image
%         scatter(tempx(1:5:end),tempy(1:5:end),0.1,'MarkerEdgeColor','none',...
%             'MarkerFaceColor',Cmap(f,:),'MarkerFaceAlpha',0.1);
        
        % calculate 2D smooth histogram
        [ctr1,ctr2,Counts] = smoothhist2D([tempx tempy],smooth_val,[nbins nbins],0,'none'); 
%        
        [Nrow,Ncol,~] = size(Counts); % number of rows and columns of histogram       
        % create a color mask same size as histogram
        Color = cat(3,repmat(simple_maps(f,1),Nrow,Ncol),repmat(simple_maps(f,2),Nrow,Ncol),...
            repmat(simple_maps(f,3),Nrow,Ncol));
        % draw the color mask
        hh = imagesc(ctr1,ctr2,Color);colormap(CmapRed);
        % change the transperancy of color mask according to 2d hist such
        % that we will see colors only when heat map is high:
        set(hh,'AlphaData',Counts);
        
        title([Titles_dark{d},' fish: ',num2str(f)]);
        % draw only a heat map next to the image for better view
        subplot(Nf,2,f*2); 
        hh = imagesc(ctr1,ctr2,Counts);colormap(CmapRed); hold on;
        set(gca,'Xtick',[],'XtickLabel',[],'Ytick',[],'YtickLabel',[]);
        axis image;
        if ~isempty(img)
            axis([0 size(img,2) 0 size(img,1)]); 
        end
        
    end
    if SAVE
       name = ['Heatmaps_',Titles_dark{d}];
       saveas(gcf,[save_folder,name]);
       disp(['saved heatmap for ',Titles_dark{d}]);
    end
end






%%  plot examples for active and inactive events:
le = 2; % light experiment
de = 1; % dark experiment
f_num = 1;  % fish number


f = figure();
subplot(1,2,1);
plot(Scl{le}(:,f_num)); hold on
backgroundBoxes(gcf,[lims_l{le}{f_num}(1,1:100); lims_l{le}{f_num}(2,1:100)],Cmap(1,:),0.3);
xlim([0 1100]); ylim([0 12])
box off; title('Light');
xlabel('Time [frames]'); ylabel('Speed [pix/frames]');

subplot(1,2,2);
plot(Scd{de}(:,f_num),'Color',Cmap(2,:)); hold on
backgroundBoxes(gcf,[lims_d{de}{f_num}(1,1:100); lims_d{de}{f_num}(2,1:100)],Cmap(2,:),0.3);
xlim([0 1100]); ylim([0 12])
box off; title('Dark');
xlabel('Time [frames]'); ylabel('Speed [pix/frames]');

if SAVE
    name = ['example_active_inactive'];
    saveas(gcf,[save_folder,name]);
    disp(['saved example for active and inactive states']);
end


%% make example movie to compare fish speeds with tracking

le = 2; % light experiment
de = 2; % dark experiment
f = 1;  % fish number

% create movie objects
movl = vid_l{le}; 
movd = vid_d{de};
movl.CurrentTime = 0;
movd.CurrentTime = 0;

if MAKEMOVIE
        movie_name = ['compare_light_dark_and_speed_fish1_smoothed'];
        obj = VideoWriter([movie_name],'MPEG-4');
        obj.FrameRate = 25;
        obj.Quality = 100;
        open(obj);
end

figure('Units','Normalized','Position',[0.0 0.0 1 1]);
pause(0.5);

Tail = 50;


for i = 1:10000
    Gl = readFrame(movl);
    Gd = readFrame(movd);
    subplot(2,2,1);
    imshow(Gl); hold on;
    plot(xcl{le}(i,f),ycl{le}(i,f),'k*','MarkerSize',10);
   title([Titles_light{le},' ',num2str(i),' Q: ',num2str(qhl{le}(i,f))])
   hold off

   subplot(2,2,3);
   if i>Tail
    plot(0-Tail:0+Tail,Shl{le}(i-Tail:i+Tail,f)); ylim([0 10]); xlim([0-Tail 0+Tail]);hold on;
    plot(0-Tail:0+Tail,Scl{le}(i-Tail:i+Tail,f),'k');
    hold off;
    box off;
    xlabel('Time from current frame [frames]');
    ylabel('Speed [pix/frame]');
   end
    subplot(2,2,2);
    imshow(Gd);hold on;
    plot(xhd{de}(i,f),yhd{de}(i,f),'k*','MarkerSize',10);
   title([Titles_dark{de},' ',num2str(i)])
    hold off
       subplot(2,2,4);
   if i>Tail
        plot(0-Tail:0+Tail,Shd{de}(i-Tail:i+Tail,f)); ylim([0 10]); xlim([0-Tail 0+Tail]);
        hold on;
        plot(0-Tail:0+Tail,Scd{de}(i-Tail:i+Tail,f),'k');
        hold off;
        box off;
        xlabel('Time from current frame  [frames]');
        ylabel('Speed [pix/frame]');
   end
   pause(0.01);
   if MAKEMOVIE && i>Tail
    frame = getframe(gcf);
    writeVideo(obj,frame);
   end
end

if MAKEMOVIE
   close(obj); 
end




%%

% 
% Meds = zeros(2,Nf);
% for j= 1:Nf
%     Meds(1,j) = median(Shl(qsl(:,j),j));
%     Meds(2,j) = median(Shd(qsd(:,j),j));
% end
% 
% figure();
% temp = Shl(:,1); temp(qsl(:,1)) = nan;
% plot(temp);
% 
% figure();
% temp = Shd(:,1); temp(qsd(:,1)) = nan;
% plot(temp);
% 
% histogram(log(Shl));
