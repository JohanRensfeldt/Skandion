function data = importfile_sigma_comparision2(workbookFile, sheetName, numberOfRows, columns_to_read, startcol)
%IMPORTFILE Import data from a spreadsheet

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    columns_to_read = [1, Inf];
end

% Calculate the number of columns to read
numberOfColumns = columns_to_read(2) - columns_to_read(1) + 1;

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", numberOfColumns);

% Specify sheet and range
startColumn = startcol;  
endColumn = startColumn + numberOfColumns - 1;  % Calculate end column
startColumnLetter = numberToExcelColumn(startColumn);  % Convert column number to letter
endColumnLetter = numberToExcelColumn(endColumn);  % Convert column number to letter
opts.DataRange = strcat(startColumnLetter, num2str(columns_to_read(1)), ":", endColumnLetter, num2str(numberOfRows));  % Define range

% Import the data
opts.Sheet = sheetName;
data = readtable(workbookFile, opts, "UseExcel", false);

%% Convert to output type
data = table2cell(data);
numIdx = cellfun(@(x) ~isnan(str2double(x)), data);
data(numIdx) = cellfun(@(x) {str2double(x)}, data(numIdx));

end

function excel_column = numberToExcelColumn(n)
    % This function converts a column number to an Excel column letter

    % Define the letters used in Excel column naming
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    num_letters = length(letters);

    % Initialize the output
    excel_column = "";

    while n > 0
        if mod(n,num_letters) == 0
            excel_column = strcat(letters(num_letters), excel_column);
            n = n / num_letters - 1;
        else
            excel_column = strcat(letters(mod(n,num_letters)), excel_column);
            n = floor(n / num_letters);
        end
    end
end


