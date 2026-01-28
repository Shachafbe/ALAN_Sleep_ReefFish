function [Sleep_pct,Total_sleep_pct] = measureSleepMinute(data,time_corrected,varargin)
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input: 



%Output: 

default_speed_th  = 1; % below this value we consider immobile
default_sleep_th  = 6; % above this value we consider sleeping
default_fps  = 12; % frames per second



% parse
vars = inputParser;
addParameter(vars,'speed_th',default_speed_th);
addParameter(vars,'sleep_th',default_sleep_th);
addParameter(vars,'fps',default_fps);

parse(vars,varargin{:})

speed_th = vars.Results.speed_th;
sleep_th = vars.Results.sleep_th;
fps = vars.Results.fps;


%...........Local Variable definitions..........

[T,N] = size(data);


T_min = (1:T)./fps/60;

max_min = floor(T_min(end)); % heighest full min in experiment
%.................Main Function.................

% for every fish mesasure to %of time sleeping out of an hour (with sloop
% defines as not moving for more then th seconds
Sleep_pct = zeros(max_min,N); 
Total_sleep_pct = zeros(N,1);
for f = 1:N
   temp = data(:,f) < speed_th;
   if iscell(time_corrected) 
       temp_good_data = ~time_corrected{f};
   else
       temp_good_data = ~time_corrected(:,f);
   end
   
   temp(~temp_good_data) = false;
   [st,nd] = findNonZeroSeq(temp);
   L = nd-st+1;
   i_sleep = L/fps > sleep_th;

   Total_sleep_pct(f) = sum(L(i_sleep))./sum(temp_good_data);
   
   
   for h = 1:max_min
       temp_h = temp(T_min > h-1 & T_min < h);
       temp_good_data_h = temp_good_data(T_min > h-1 & T_min < h);
       temp_h(~temp_good_data_h) = false;
       [st,nd] = findNonZeroSeq(temp_h);
       L = nd-st+1;
       % sleep events:
       i_sleep = L/fps > sleep_th;
       Sleep_pct(h,f) = sum(L(i_sleep))./sum(temp_good_data_h);
       
   end

end








%............Call for local functions...........




