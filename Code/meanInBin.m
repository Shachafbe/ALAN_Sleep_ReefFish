function [M,binc,nn,S,l2] = meanInBin(data,responses,bins,varargin)
%description: This function gets a continuous variables data, and a
%continuous response vector of the same length, and calculates the average
%response values (and std) in each bin. Can also return median and iqr.
% Predictor data can be either 1d or 2d. 

%Input: data - nx1 or nx2 matrix of predictor data. responses nx1 vecotr.
        % Optional inputs: Type - mean (std) or median (iqr)

%Output:


% set defaults
defaultPLOT = 1;
defaultType = 'mean';
defaultColorInd = 1;
defaultLabels = {'x','y'};
defaultSmooth = [];
defaultNo_data_val = 0;
defaultTH_per_bin = [];
defaultSig_per_bin = [];
defaultCmap = myCmap1; 

% parse
vars = inputParser;
addParameter(vars,'Plot',defaultPLOT);
addParameter(vars,'Labels',defaultLabels);
addParameter(vars,'Smooth',defaultSmooth);
addParameter(vars,'No_data_val',defaultNo_data_val);
addParameter(vars,'TH_per_bin',defaultTH_per_bin);
addParameter(vars,'Sig_per_bin',defaultSig_per_bin);
addParameter(vars,'ColorInd',defaultColorInd);
addParameter(vars,'Type',defaultType);
addParameter(vars,'Cmap',defaultCmap);


parse(vars,varargin{:})

% make variables
PLOT = vars.Results.Plot;
Labels = vars.Results.Labels;
Smooth = vars.Results.Smooth;
No_data_val = vars.Results.No_data_val;
TH_per_bin = vars.Results.TH_per_bin;
Sig_per_bin = vars.Results.Sig_per_bin;
ColorInd = vars.Results.ColorInd;
Type = vars.Results.Type;
Cmap = vars.Results.Cmap;



%...........Local Variable definitions..........









%.................Main Function.................

% one dimensional data set:

if size(data,2)==1
    
    [nn,bins,ii] = histcounts(data,bins); % bin the data
    bin_width = (bins(2)-bins(1)); % calculate bin width
    binc = bins(1:end-1)+bin_width/2; % calc bin centers
    M = zeros(1,length(nn)); % calcualte averages
    S = M;
    for j = 1:length(M)
        temp = find(ii==j);
        if ~isempty(temp)
            if strcmp(Type,'mean')
                M(j) = mean(responses(temp));
                S(j) = std(responses(temp));
            elseif strcmp(Type,'median')
                M(j) = nanmedian(responses(temp));
                S(j) = iqr(responses(temp))/2;
            elseif strcmp(Type,'sum')
                M(j) = sum(responses(temp));
                S(j) = 0;
            else
                disp('type doesnt match')
                return
            end
                
        end
    end
    

    
    % plot if needed:
    if PLOT
        l2 = boundedLine(binc,M,S);
        if exist('Labels','var')
            xlabel(Labels);
        end
    end
    
    % if we have two dimensional data set:
elseif size(data,2)==2
    
    % bin the fisrt variable
    [N1,E1] = discretize(data(:,1),bins{1});
    
    % bin the second variable
    [N2,E2] = discretize(data(:,2),bins{2});
    
    % calculate probability to bout in each bin
    M = zeros(max(N1),max(N2));
    S = M;
    nn = S;
    for i = 1:size(S,1)
        for j = 1:size(S,2)
            temp = N1==i & N2==j;
            nn(i,j) = sum(temp); % save the number of jount events
            % and the probability to bout:
            if sum(temp)>0
                if strcmp(Type,'mean')
                    M(i,j) = mean(responses(temp));
                    S(i,j) = std(responses(temp));
                elseif strcmp(Type,'median')
                    M(i,j) = median(responses(temp));
                    S(i,j) = iqr(responses(temp));
                elseif strcmp(Type,'sum')
                    M(i,j) = sum(responses(temp));
                    S(i,j) =0;
                else
                     disp('type doesnt match')
                     return
                end
            end
        end
    end
    

    if PLOT
        figure();
        M(nn<TH_per_bin) = 0;
        S(nn<TH_per_bin) = 0;
        subplot(2,2,1);
        l2(1) = imageSC(E1,E2,M');colorbar; colormap(CmapRed); axis square;
        title('Mean property');
        if exist('Labels','var')
            xlabel(Labels{1});ylabel(Labels{2});
        end
        Lims = evenLimits(M(~isinf(M)));
        colormap(Cmap);caxis(Lims);
        subplot(2,2,2);
        l2(2) = imageSC(E1,E2,nn'); colorbar; axis square;
        if exist('Labels','var')
            xlabel(Labels{1});ylabel(Labels{2});
        end
        title('Bin count');
        subplot(2,2,3);
        l2(3) = imageSC(E1,E2,S');colorbar; axis square;
        if exist('Labels','var')
            xlabel(Labels{1});ylabel(Labels{2});
        end
        title('Error size');
        caxis(evenLimits(S(~isinf(S))));
    end
end
%............Call for local functions...........




