function [Sleep_pct,Total_sleep_pct] = measureSleep(data,time_corrected,varargin)
%description: 


%Input: 



%Output: 

default_speed_th  = 1; % below this value we consider immobile
default_sleep_th  = 6; % above this value we consider sleeping
default_sleep_per_time  = 'hour'; % 'hour' or 'min'.
default_fps  = 22; % frames per second



% parse
vars = inputParser;
addParameter(vars,'speed_th',default_speed_th);
addParameter(vars,'sleep_th',default_sleep_th);

addParameter(vars,'sleep_per_time',default_sleep_per_time);

addParameter(vars,'fps',default_fps);

parse(vars,varargin{:})

speed_th = vars.Results.speed_th;
sleep_th = vars.Results.sleep_th;
fps = vars.Results.fps;

sleep_per_time = vars.Results.sleep_per_time;

%...........Local Variable definitions..........

[T,N] = size(data);

if strcmp(sleep_per_time,'hour') % if it's sleep per hour
    T_hr = (1:T)./fps/60/60;
elseif strcmp(sleep_per_time,'min') % otherwise it's per minute
    T_hr = (1:T)./fps/60;
else
    disp('Error: sleep_per_time should be either ''hour'' or ''min'' ');
    return
end

max_hr = floor(T_hr(end)); % heighest full hour in experiment
%.................Main Function.................

% for every fish mesasure to %of time sleeping out of an hour (with sloop
% defines as not moving for more then th seconds
Sleep_pct = zeros(max_hr,N); 
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
   
   
   for h = 1:max_hr
       temp_h = temp(T_hr > h-1 & T_hr < h);
       temp_good_data_h = temp_good_data(T_hr > h-1 & T_hr < h);
       temp_h(~temp_good_data_h) = false;
       [st,nd] = findNonZeroSeq(temp_h);
       L = nd-st+1;
       % sleep events:
       i_sleep = L/fps > sleep_th;
       Sleep_pct(h,f) = sum(L(i_sleep))./sum(temp_good_data_h);
       
   end

end








%............Call for local functions...........




