function [output] = watchSingleFishMovie(moviename,Start,End,x,y,varargin)
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input: 



%Output: 

default_var  = [];



% parse
vars = inputParser;
addParameter(vars,'exmp_var',default_var);

parse(vars,varargin{:})

exmp_var = vars.Results.exmp_var;


%...........Local Variable definitions..........

vid = VideoReader(moviename);
Fs = vid.FrameRate;


%.................Main Function.................
figure('Units','Normalized','Position',[0.01 0.2 0.8 0.7]);
for i = Start:End
    vid.CurrentTime = (i/Fs);
    fr = readFrame(vid);
    imagesc(fr); axis image; colormap gray; hold on;
    plot(x(i),y(i),'r.');
    title(['frame: ',num2str(i), 'Time: ',num2str(i/Fs,3),' [s]']);
    pause(0.01);
    hold off;
end



%............Call for local functions...........




