function f = backgroundBoxes(f,xx,Cmap,Alpha)
%This funciton adds background boxes to a figure

%Input: xx - a matrix of 2xp entries (2-is start and end of each box p-number of polygons)
%       f - is a handle for a figure
%       Cmap - a matrix of Nx3 - which is the number of different colors to
%              use
%       Alpha - sets the transperancy of the boxes (optional);

%Output:



%...........Local Variable definitions..........
if ~exist('Alpha','var')
    Alpha = 1;
end

h = f.CurrentAxes;
%get number of colors
N = size(Cmap,1);
p = size(xx,2);

%get Y upper limit and lower limit
Ymax = h.YLim(2);
Ymin = min([h.YLim(1) 0]);
yy = [Ymin Ymax Ymax Ymin];

% set color indices
cind = mod(1:p,N); cind(cind==0) = N;
%.................Main Function.................


%reshape the data to get boxes (original data gives the x-positions, of start)
figure(f);
%plot boxes to current figure handle
for i = 1:p
    tempx = [xx(1,i) xx(1,i) xx(2,i) xx(2,i)];
    pp = patch(tempx,yy,Cmap(cind(i),:),'EdgeColor','none');
    if Alpha<1
        pp.FaceAlpha = Alpha;
    end
end



%............Call for local functions...........




