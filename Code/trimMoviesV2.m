function trimMovies
%description: 

%Input: 



%Output: 




%...........Local Variable definitions..........

% path to file location .xlsx
Path_location_xlsx = 'G:\My Drive\PhD\Stimuli Data\RAW DATA old tracking\table.xlsx'; %

location_file = readtable(Path_location_xlsx);

L = height(location_file); % length of files:

% path to stimulus times files
Path_stim_times_xlsx ='G:\My Drive\PhD\Stimuli Data\RAW DATA old tracking\fish_stim_times.csv';

stim_times_file = readtable(Path_stim_times_xlsx);


Win_size = 90; % time in sec before and after stimuli



% if ~exist(savepath,'dir')
%     mkdir(savepath);
%     disp([savepath ,' created']);
% else
%     disp([savepath ,' exists']);
% end



for i = 1:L % loop over all files
    try
    curr_file = location_file.File(i); curr_file = curr_file{1};
    Arena = location_file.Arena(i); Arena = str2num(Arena{1}(end));
    curr_lux = location_file.LUX(i); 
    curr_fish = location_file.Fish(i); 
    
    curr_fish = regexp(curr_fish{1},'\d*','Match'); 
    
    if ~isempty(curr_fish)
        curr_fish = str2double(curr_fish{end});
    end
    curr_trial = location_file.Trial(i);
    curr_trial = regexp(curr_trial{1},'\d*','Match'); 
    curr_trial = str2double(curr_trial{end});
    % if these are tracked fish
    if Arena == 1 && ~isempty(curr_fish)
       % find start time in file:
       curr_trial_i = stim_times_file.fishnumbers==curr_fish & stim_times_file.trialnumbers==curr_trial;
       if isempty(curr_trial_i) || curr_trial_i==0
           disp(['cant find fish stim time, skipping file ',curr_file])
           continue
       end
       curr_time = stim_times_file.stimtimes{curr_trial_i};
       
       if ~isempty(curr_time)
           % transform time to seconds:
           F = 'mm:ss';
           du = duration(curr_time, 'InputFormat', F);
           stim_time_in_sec= seconds(du);
           
       else
          disp(['no time for file ',curr_file]);
          continue
       end
       
       % read movie file
       obj = VideoReader(curr_file);
       fps = obj.FrameRate;
        
       % create new movie
       movie_name = [curr_file(1:end-4),'_short'];
       vid = VideoWriter(movie_name,'Motion JPEG AVI');
       vid.FrameRate = fps;
       open(vid);
    
        curr_win = min([stim_time_in_sec,Win_size]);
        num_frames_in_movie = curr_win*2*fps;
        curr_time = stim_time_in_sec - curr_win;
        obj.CurrentTime = curr_time;
    
        cnt = 1;
        while cnt <= num_frames_in_movie
            fr = readFrame(obj);
            %         fr = rgb2gray(fr);
            writeVideo(vid,fr);
            cnt = cnt+1;
        end
        
        disp(['finished file ',curr_file]);
        
        close(vid)
        
    end
    catch 
        disp(['error for file ',curr_file])
    end
end
    
    
   
    
    
   
    
    
    


    
    


    
    
    
    
    
    
    


%.................Main Function.................



%............Call for local functions...........




