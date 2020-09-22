function [edat] = import_ecopuck_raw(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [DATE1,TIMER,BB532,BB700,CHL,OTHER] = IMPORTFILE(FILENAME) Reads data
%   from text file FILENAME for the default selection.
%
%   [DATE1,TIMER,BB532,BB700,CHL,OTHER] = IMPORTFILE(FILENAME, STARTROW,
%   ENDROW) Reads data from rows STARTROW through ENDROW of text file
%   FILENAME.
%
% Example:
%   [date1,timer,bb532,bb700,chl,other] = importfile('EcoPuck_CTD001_201805242122.raw',2, 6219);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2018/05/24 22:42:14

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf-1;
end

%% Format for each line of text:
%   column1: text (%q)
%	column2: text (%q)
%   column4: double (%f)
%	column6: double (%f)
%   column8: double (%f)
%	column9: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%*q%f%*q%f%*q%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
date1 = dataArray{:, 1};
timer = dataArray{:, 2};
edat.bbs532.raw = dataArray{:, 3};
edat.bbs700.raw = dataArray{:, 4};
edat.chl.raw = dataArray{:, 5};
edat.other = dataArray{:, 6};

%% calculate concentrations
edat.bbs532.scale_factor=6.974*10^-6;
edat.bbs532.dark_counts=53;
edat.bbs532.concentration=edat.bbs532.scale_factor*(edat.bbs532.raw-edat.bbs532.dark_counts);

edat.bbs700.scale_factor=3.004*10^-6;
edat.bbs700.dark_counts=52;
edat.bbs700.concentration=edat.bbs700.scale_factor*(edat.bbs700.raw-edat.bbs700.dark_counts);

edat.chl.scale_factor=.0305;
edat.chl.dark_counts=51.5;
edat.chl.concentration=edat.chl.scale_factor*(edat.chl.raw-edat.chl.dark_counts);

aa=nansum(~ismissing(char(timer)));
for i=1:aa(1)
    d=char(date1(i));
    t=char(timer(i));
    time(i)=datenum(str2double(d(7:8))+2000,str2double(d(1:2)),str2double(d(4:5)),str2double(t(1:2)),str2double(t(4:5)),str2double(t(7:8)));
end

edat.time=time';