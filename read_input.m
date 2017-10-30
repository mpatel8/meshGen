% clc; clear all ;

function [mN, mE] = read_input (saFlInput)
%saFlInput = 'seg-1-36.inp' ;
fid = fopen(saFlInput,'rt') ;
TXT = textscan(fid,'%s','Delimiter','\n');
TXT = TXT{1} ;
%% Get the line number of mises 
idxN = strfind(TXT, 'Node');
idx1 = find(not(cellfun('isempty', idxN)));
idxE = strfind(TXT, 'Element');
idx2 = find(not(cellfun('isempty', idxE)));
if isempty(idx2)
    idxE = strfind(TXT, 'ELEMENT');
    idx2 = find(not(cellfun('isempty', idxE)));
end
idx2 = find(not(cellfun('isempty', idxE)));
idxS = strfind(TXT, 'End');
idx3 = find(not(cellfun('isempty', idxS)));
% pick  nodes 
nodes = TXT(idx1+1:idx2-1);
mN = cell2mat(cellfun(@str2num,nodes,'UniformOutput',false));

% pick elements 
elements = TXT(idx2+1:idx3(1)-1);
count = 0;
ele = cell(length(elements)/2,1);
mE = cell2mat(cellfun(@str2num,elements,'UniformOutput',false));
% % for i = 1:2:length(elements)
% %     count = count+1;
% %     ele{count,1} = [elements{i} elements{i+1}];
% % end
% % mE = cell2mat(cellfun(@str2num,ele,'UniformOutput',false));
