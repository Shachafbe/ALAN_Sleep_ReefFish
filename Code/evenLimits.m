function Limits = evenLimits(X,Reminder)
%function gets a vector or a matrix of data and calculate the maximal and
%minimum values so that zero will be the middle value.

%Input: X - matrix or vector

%Output: Limits - lower and upper limits
if nargin==1
    Reminder = 0;
end

X = X(:); %make a col vec

%get the maximum abolute difference from zero;
temp = max(abs([max(X) min(X)]));

% add some reminder to the limits in a form of added percentage (Default is 0)
Addition = Reminder*temp;
temp = temp+Addition;

%set it as limits
Limits = [-temp temp];



%...........Local Variable definitions..........



%.................Main Function.................



%............Call for local functions...........




