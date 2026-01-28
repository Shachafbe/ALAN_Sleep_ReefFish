function [output] = runExtractSingleFish
%description: allows user to draw ROI over a brain and extract the
%identity of cells included in that ROI


%Input: 



%Output: 
FileNames = {'314.avi','314.avi'};
FishNames = {'fish7_314','fish8_314'};


for i = 1:length(FileNames)
    extractObjectsSingleFishMovieFiles(FileNames{i}, FishNames{i},1);
    disp('finished ',FishNames{i})
end


%...........Local Variable definitions..........


%.................Main Function.................



%............Call for local functions...........




