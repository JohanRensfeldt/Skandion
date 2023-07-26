function GTR1andGTR2SpotSizeG315z0sigmacomparision = importfile_sigma_comparision(workbookFile, sheetName, dataLines)
%IMPORTFILE Import data from a spreadsheet
%  GTR1ANDGTR2SPOTSIZEG315Z0SIGMACOMPARISION = IMPORTFILE(FILE) reads
%  data from the first worksheet in the Microsoft Excel spreadsheet file
%  named FILE.  Returns the data as a cell array.
%
%  GTR1ANDGTR2SPOTSIZEG315Z0SIGMACOMPARISION = IMPORTFILE(FILE, SHEET)
%  reads from the specified worksheet.
%
%  GTR1ANDGTR2SPOTSIZEG315Z0SIGMACOMPARISION = IMPORTFILE(FILE, SHEET,
%  DATALINES) reads from the specified worksheet for the specified row
%  interval(s). Specify DATALINES as a positive scalar integer or a
%  N-by-2 array of positive scalar integers for dis-contiguous row
%  intervals.
%
%  Example:
%  GTR1andGTR2SpotSizeG315z0sigmacomparision = importfile("C:\Users\G_rej001\Desktop\Plot_Data\GTR1 and GTR2 SpotSize_G315_z0_sigma_comparision.xlsx", "GTR1 G315 z0 x", [9, 52]);

%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 2
    dataLines = [9, 52];
end

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 53);

% Specify sheet and range
opts.Sheet = sheetName;
opts.DataRange = "D" + dataLines(1, 1) + ":BD" + dataLines(1, 2);

% Specify column names and types
opts.VariableNames = ["VarName4", "REF", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56"];
opts.VariableTypes = ["char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char"];

% Specify variable properties
opts = setvaropts(opts, ["VarName4", "REF", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName4", "REF", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55", "VarName56"], "EmptyFieldRule", "auto");

% Import the data
GTR1andGTR2SpotSizeG315z0sigmacomparision = readtable(workbookFile, opts, "UseExcel", false);

for idx = 2:size(dataLines, 1)
    opts.DataRange = "D" + dataLines(idx, 1) + ":BD" + dataLines(idx, 2);
    tb = readtable(workbookFile, opts, "UseExcel", false);
    GTR1andGTR2SpotSizeG315z0sigmacomparision = [GTR1andGTR2SpotSizeG315z0sigmacomparision; tb]; %#ok<AGROW>
end

%% Convert to output type
GTR1andGTR2SpotSizeG315z0sigmacomparision = table2cell(GTR1andGTR2SpotSizeG315z0sigmacomparision);
numIdx = cellfun(@(x) ~isnan(str2double(x)), GTR1andGTR2SpotSizeG315z0sigmacomparision);
GTR1andGTR2SpotSizeG315z0sigmacomparision(numIdx) = cellfun(@(x) {str2double(x)}, GTR1andGTR2SpotSizeG315z0sigmacomparision(numIdx));
end