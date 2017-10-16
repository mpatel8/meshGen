% clc; clear all ;

function [mN, mE] = read_input (saFlInput)
%saFlInput = 'seg-1-36.inp' ;
fid = fopen(saFlInput,'rt') ;
TXT = textscan(fid,'%s','Delimiter','\n');
TXT = TXT{1} ;
%% Get the line number of mises 
idxS = strfind(TXT, 'Node');
idx1 = find(not(cellfun('isempty', idxS)));
idxS = strfind(TXT, 'Element');
idx2 = find(not(cellfun('isempty', idxS)));
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
