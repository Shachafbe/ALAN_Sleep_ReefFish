function [binned_prop,cont_binned_prop, tt] = plotContProperty(prop,time_bin,varargin)
%description: bins the data stored in prop, rarranges it in continious time, and plots 


%Input: prop - cell array with property to plot
%       Lframes - number of frames in each cell of prop
%       time_bin - time in minutes to bin the data
%       Fs - frame rate



%Output: % binnd_prop 
%        % cont_binned_prop - continious according to [dim1,(dim3,dim2,dim4)]



%...........Local Variable definitions..........
% set defaults
default_Fs =  1; %frame rate
default_prop_i = 1; % property idx to plot
defaultPLOT = 0;
defaultSAVE = 0;




% parse
vars = inputParser;
addParameter(vars,'Fs',default_Fs);
addParameter(vars,'prop_i',default_prop_i);
addParameter(vars,'PLOT',defaultPLOT);
addParameter(vars,'SAVE',defaultSAVE);



parse(vars,varargin{:})
% make variables
Fs =  vars.Results.Fs; %m
prop_i =  vars.Results.prop_i; %m
PLOT = vars.Results.PLOT;
SAVE = vars.Results.SAVE;


%.................Main Function.................
%% plot speed over time for all fish 
Lframes = cellfun(@length,prop);

temp_Lframes = Lframes(:);

min_frames = min(temp_Lframes(temp_Lframes~=0));
% time_bin = 20; % in min
bin_size = time_bin*(60*Fs); % bin size in minutes
time_bins = 0:bin_size:min_frames;

% varialbe for binned data
binned_prop = zeros([size(prop) length(time_bins)-1]);

for f = 1:size(prop,1)
    for i = 1:size(prop,2)
        for j = 1:size(prop,3)
            for m = 1:size(prop,4)
                if ~isempty( prop{f,i,j,m})
                temp_s = prop{f,i,j,m}(:,prop_i);
                
                % bin the data 
                [M,binc,nn,S] = meanInBin((1:length(temp_s))',temp_s,time_bins,'PLOT',0);
                binned_prop(f,i,j,m,:) = M;
                end
            end
        end
    end
end

% now plot the data over all fish

% first rearrange as continious
cont_binned_prop = zeros(size(prop,1),(size(binned_prop,5)-1)*8)*nan;
for f = 1:size(binned_prop,1)
    p = 1;
    for c = 1:size(binned_prop,2)
        for d = 1:size(binned_prop,3)
            for t = 1:size(binned_prop,4)
                
                temp = squeeze(binned_prop(f,d,c,t,:));
                cont_binned_prop(f,p:p+length(temp)-1) = temp;
                p = p+length(temp);
                
            end
        end
    end
end

% cont_binned_prop(cont_binned_prop==0) = nan;
tt = linspace(1,12*8,length(cont_binned_prop)); % time in hours

if PLOT
figure();
plot(tt,cont_binned_prop','LineWidth',1,'color',[0 0 0 0.5]);
M = nanmean(cont_binned_prop);
stde = nanstd(cont_binned_prop)./sqrt(sum(~isnan(cont_binned_prop)));
boundedLine(tt,M,stde);
st = linspace(0,max(tt),9);
nd = st(2:end);
st(end) = [];
xx = [st;nd];
cmap1 = [1 1 1; 0 0 0];
cmap2 = [0.7 0.7 0.7; 0 0 0 ];
backgroundBoxes(gcf,xx(:,1:4),cmap1,0.2);
backgroundBoxes(gcf,xx(:,5:8),cmap2,0.2);
box off; 

xlim([0 12*8])
end

%............Call for local functions...........




