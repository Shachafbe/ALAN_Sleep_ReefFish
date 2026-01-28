function [Strt,End, vec] = findNonZeroSeq(vec,LowerLimit,connect_th)
%this function gets a vectos of scalars and returns beginning and end positions of
%non zero sequences in it

% optional input LowerLimts - if given - removes all events shorter or
% equal to Lower limit and return a vector that replace these events with zeros 

%make sure we have a row vector
isrow = 1;
if size(vec,1) > size(vec,2)
    isrow = 0;
    vec = vec';
end

% vec = [0 1 1 1 0 0 1 1 1 1];
limits = diff([0, vec, 0]);
Strt = find(limits==1);
End = find(limits==-1)-1;


% remove short events if needed:
if exist('LowerLimit','var') && ~isempty(LowerLimit)
    L = End-Strt + 1;
    
    % remove short events
    ii_remove = find(L<LowerLimit);
    for i = 1:length(ii_remove)
        vec(Strt(ii_remove(i)):End(ii_remove(i))) = 0;
    end
    
    limits = diff([0, vec, 0]);
    Strt = find(limits==1);
    End = find(limits==-1)-1;
    
    if ~isrow
        vec = vec';
    end
    
end



% connect events if needed
if exist('connect_th','var')
    
    while connect_th > 0
    
    
        % find distances between event
        Lbetween = Strt(2:end)-End(1:end-1);
        
        i_to_connect = find(Lbetween<connect_th);
        
        if isempty(i_to_connect)
            connect_th = 0;
        end
        
        for i = 1:length(i_to_connect)
            
            curr_i = i_to_connect(i);
            
            vec(End(curr_i):Strt(curr_i+1)) = 1;
            
        end
            End(i_to_connect) = [];
            Strt(i_to_connect+1) = [];
    end
    if ~isrow
        vec = vec';
    end
    
end

