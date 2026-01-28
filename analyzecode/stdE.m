function stde = stdE(data,dim)
%This function calculates the std error of a data set X
%If "data" is a matrix it calculates the std error of each column.
%Input: data - a matrix(or vector) of data
%Output: stde - a vector containing the std error of the rows or columns
%depending on dim value

if nargin<2
    dim = 1; %default is according to row
end

% if size(data,1)==1 %if data is a row vector make a cloumn
%     data = data';
% end

%calculate std error according to dim given
if sum(isnan(data(:))) == 0
    
    st_d = std(data,[],dim);
    %devide std by the squre root of the length of the vector
    stde = st_d/sqrt(size(data,dim));
    
else % add calculation in case of nans
    st_d = nanstd(data,[],dim);
    N = sum(~isnan(data),dim);
    
    stde = st_d./sqrt(N);
    
end
