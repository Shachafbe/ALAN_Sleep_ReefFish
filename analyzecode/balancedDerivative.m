function Deriv = balancedDerivative(data,varargin)
%description: calculates a balanced derivative such that dx(t) =
% x(t+L)-x(t-L)/2L and padds with zeros s.t. the length of the outputs is
% the same as the input


%Input: data array if derivative is 1d or x,y id derivative is 2d
%       optional arguments:
%       dim - dimensiotn to calculate derivative (1,2 or 3)
%       L - length of derivative (2L)


%Output: Deriv - same size as data



%...........Local Variable definitions..........

% set defaults
defaultDim = 1;
defaultL = 1;

% parse
vars = inputParser;
addParameter(vars,'dim',defaultDim);
addParameter(vars,'L',defaultL);


parse(vars,varargin{:})
% make variables
dim = vars.Results.dim;
L = vars.Results.L;


%Input:


%.................Main Function.................
if dim==1
    % take difference with distance 2*L
    temp = (data(L+2:end,:) - data(1:end-(L+1),:))/(2*L);
    % padd with zeros
    Deriv = [temp(1,:); temp; temp(end,:)];
    
    
elseif dim == 2
     % take difference with distance 2*L
    temp = (data(:,L+2:end) - data(:,1:end-(L+1)))/(2*L);
    % padd with zeros
    Deriv = [temp(:,1) temp temp(:,end)];
    
elseif dim == 3
    temp = (data(:,:,L+2:end) - data(:,:,1:end-(L+1)))/(2*L);

    Deriv = cat(3,temp(:,:,1),temp,temp(:,:,end));

end

%............Call for local functions...........




