clearvars
close all 

choice = questdlg('Would you like to include outliers?', 'Include Outliers', 'Yes', 'No', 'No');
include_outliers = strcmp(choice, 'Yes');

choice = questdlg('Would you like to plot several energies?', 'Plot several', 'Yes', 'No', 'No');
plot_several = strcmp(choice, 'Yes');

%% Calling functions and importing data
% This section determines the number of columns in the data file and
% then imports the data. Here the user also chooses what dates to look at 
% and what energies.


[file, path] = uigetfile('*.xlsx','Select an Excel data file');
if isequal(file, 0)
    disp('User selected Cancel')
else
    disp(['User selected ', fullfile(path, file)])
    % Get all sheet names
    allSheets = sheetnames(fullfile(path, file));
    % Let the user select a sheet
    [indx,tf] = listdlg('ListString', allSheets, 'SelectionMode', 'single');
    if tf == 0
        disp('User selected Cancel')

    else
        disp(['User selected ', allSheets{indx}])
        sheet = allSheets{indx};
        path = fullfile(path, file);
        % Load data from the selected sheet
        number_cols = number_of_columns(path, sheet);
        data = importfile_sigma_comparision2(path, allSheets{indx}, 52, [1, number_cols], 4);  % Modify this as needed to match your data format       
    end
end

data = data(9:29, :);

start_date = '10-Jan-2019' ;

end_date = '20-Jun-2023' ;

if strcmp(sheet, 'GTR2 G315 z0 x') || strcmp(sheet, 'GTR2 G315 z0 y')

    start_date = '22-Oct-2019' ;

    end_date = '27-May-2023' ;

end

% [start_date, end_date] = select_range() ;

mev = 60 ;

mev_several = [60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 226];

[start_col, end_col] = match_date_to_column(data,start_date,end_date);

datetimes = convertDates(data, start_col, end_col);

data_Table = covert_to_data_table(data) ;

table_averages = yearly_averages(data_Table); 

make_boxplot(data_Table, sheet)

if plot_several
    diffrences_array = plot_several_energies(data,start_col, end_col, mev_several, datetimes, include_outliers, sheet);
else
    plot_data(data,start_col,end_col, mev, datetimes, include_outliers, sheet)
end

create_trend_lines(diffrences_array, sheet, mev_several)

plot_average_data(table_averages, sheet)

%% Plotting the averages of the data
% In this function the script takes in the avrages of the energy levels
% over the data range. The script then calculates the error relative to the
% reference value and then fits a second degree plynomial to the data and
% plots the results.


function plot_average_data(averages,sheet)
    
    reference = averages(:, 2);
    
    data = averages{:, 2:end};

    averagesRows = mean(data, 2);
    
    ener = [60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 226];

    relative_error = 1 - reference ./ averagesRows; 
    
    relative_error_a = table2array(relative_error) * 100;

    p = polyfit(ener, averagesRows, 2);
    
    p1 = polyfit(ener, relative_error_a, 2);

    xRange = linspace(min(ener), max(ener), 18);

    yFit = polyval(p, xRange);

    yFit1 = polyval(p1, xRange);

    figure

    scatter(ener, averagesRows, 'b', 'filled')

    hold on

    plot(ener, yFit, 'r')
    
    title(sprintf('Spot size as a function of energy in MeV for %s',sheet))
    xlabel('Energy in Mev');
    ylabel('Spot size');
    legend('Original Data', 'Fitted curve')

    str = sprintf('y = %.6f*x^2 + %.2f*x + %.2f', p(1), p(2), p(3));
    text(max(ener)/2, min(averagesRows), str, 'FontSize', 14);

    hold off

    %%%%%%%%%%%%%

    figure

    scatter(ener, relative_error_a, 'magenta', 'filled')

    hold on

    plot(ener, yFit1, 'b')
    
    xlabel('Energy in Mev');
    ylabel('Relative deviation from reference in %');
    title(sprintf('Relative spot size deviation %s', sheet));

    legend('Mean of relative deviation from reference in %', 'Fitted curve')

    str1 = sprintf('y = %.6f*x^2 + %.6f*x + %.4f', p1(1), p1(2), p1(3));
    text(max(ener)/2, min(relative_error_a), str1, 'FontSize', 14);

    hold off

end

%% Makes boxplot, yearly averages and the standard deviation of the data.
% In this function the script creates a boxplot from the average deviation
% of the data per energy level. The function also creates a plot with the
% standard deviation and the mean of the relative error.

function make_boxplot(total_data, sheet)
    % Convert the average values in the table to a matrix
    total_data_array = table2array(total_data);

    average = zeros(height(total_data),width(total_data)-3);

    for row = 1 : height(total_data)
        for columns = 3 : width(total_data)
            ref = total_data_array(row, 2);
            data = total_data_array(row, columns);
            average(row, columns - 2) = (data/ref) - 1;
        end
    end
    
    mean_value = mean(average') * 100 ;

    std_data = std(average') * 100 ;

    average = (average * 100)' ;

    % Create labels from the energy levels
    energy_levels = total_data.("E [MeV]");
    labels = (cellstr(num2str(energy_levels)))';  % Convert to cell array of character vectors
    
    % Create a boxplot of the average data
    figure()
    boxplot(average, labels);  % Use the labels in the boxplot function
    title(sprintf('Yearly Averages for %s', sheet));
    xlabel('Energy Level (E [MeV])');
    ylabel('Relative deviation from reference (%)');
    hold off
    
    figure()
    scatterHandle1 = scatter(energy_levels, mean_value,'filled','square');
    hold on
    plot(energy_levels, mean_value,'--k');
    
    scatterHandle2 = scatter(energy_levels, std_data,'filled','o','red');
    plot(energy_levels, std_data, '--r');
    ylabel('percent deviation from reference')
    xlabel('Energy Mev')
    title(sprintf('Standard deviation and mean deviation from reference in percent for %s', sheet))
    
    legend([scatterHandle1, scatterHandle2], {'Mean value', 'Standard deviation'}, 'Location','northwest');
    hold off

end

%% Creates Yearly averages
% This function creates yearly averages from the total data in the selected
% data range.


function average_table = yearly_averages(data_table_total)
    % Assuming you have a datatable with dates as column headers and energy values as rows

    % Get the date column headers
    dates = data_table_total.Properties.VariableNames(3:end);

    % Convert the dates to datetime format
    dates = datetime(dates, 'InputFormat', 'dd-MMM-yyyy');

    % Determine the start and end dates for each one-year interval
    start_date = dates(1);
    start_date = datetime(start_date.Year, 1, 1); % Start from the beginning of the first year
    end_date = dates(end);
    end_date = datetime(end_date.Year, 12, 31); % End at the end of the last year

    start_dates = start_date:years(1):end_date;
    end_dates = start_dates + years(1);

    % Determine the size of the average table
    num_intervals = length(start_dates);
    num_rows = size(data_table_total, 1);

    % Extract the first two columns
    first_two_columns = data_table_total(:, 1:2);

    % Initialize the average table with NaN values and temporary unique variable names
    tempVarNames = "Temp" + string(1:num_intervals);
    average_table = array2table(NaN(num_rows, num_intervals), 'VariableNames', tempVarNames);

    % Add the first two columns to the average table
    average_table = [first_two_columns, average_table];

    % Calculate the average energy value for each interval
    for i = 1:num_intervals
        start_date = start_dates(i);
        end_date = end_dates(i);

        % Extract the columns within the current interval
        mask = dates >= start_date & dates < end_date;
        interval_data = data_table_total(:, mask);

        % Convert the interval_data table into an array
        interval_data_array = table2array(interval_data);

        % Calculate the average energy value for the interval
        % Exclude the first two columns ('E [MeV]' and 'REF')
        interval_avg = mean(interval_data_array(:, 3:end), 2);

        % Assign the interval average values to the average table
        average_table{:, i+2} = interval_avg; % Adjusted to account for the two additional columns
    end

    % Replace the temporary variable names with the correct year names
    % Converting dates to string format to use as table variable names
    string_dates = cellstr(datestr(start_dates, 'dd-mmm-yyyy'));
    average_table.Properties.VariableNames(3:end) = string_dates;

end

%% Creates data table
% This function takes in the data as an array and converts it to a data
% table, the varnames are the dates the data was meassured.


function data_as_table = covert_to_data_table(in_data)
    % Assume 'C' is your cell array
    
    % Combine the relevant elements from the first and second row to form headers
    varNames = [in_data(1,1), in_data(2,2:end)];
    
    varNames = matlab.lang.makeUniqueStrings(varNames);

    varNames{2} = 'REF';
    
    % Remove those rows from your cell array
    in_data(1:3,:) = [];
    
    % Convert cell array to table using the variable names
    data_as_table = cell2table(in_data, 'VariableNames', varNames);
end

%% Finds start and and column
% This function takes in the data range as a start and end date and
% converts that to a range given in the a start column and a end column.


function [col_start, col_end] = match_date_to_column(data,start_date,end_date)

    match_start_date = strcmp(data(2,:),start_date) ; 

    match_end_date = strcmp(data(2,:),end_date) ;

    col_start = find(match_start_date);

    col_end = find(match_end_date);

    if isempty(col_start)
        sprintf("Date: %s specified is of invalid format", start_date)        
    elseif isempty(col_end)
        sprintf("Date: %s specified is of invalid format", end_date)
    end
end

%% Plot singular energy level
% This function plots the data of a singular energy level with and without
% outlayers.


function plot_data(data,col_start, col_end, mev, datetimes, include_outliers)

    index_mev = cellfun(@(x) isequal(x, mev), data);
    
    data_to_plot = zeros(1,col_end-col_start);
    
    double_array = cell2mat(cellfun(@double, data(index_mev,2:end), 'UniformOutput', false));

    for i = col_start : col_end - 3 
        data_to_plot(1,i-col_start+1) =  double_array(1,1) - double_array(1,i);
    end
    
    difference_procent = (data_to_plot./double_array(1,1)) * 100;
    
    if ~include_outliers
        data_table = table(datetimes', difference_procent', 'VariableNames', {'Date', 'Data'});
        
        threshold = 10;
        
        is_above_threshold = data_table(data_table.Data > threshold, :) ;
        
        figure

        scatter(is_above_threshold.Date, is_above_threshold.Data)
        
        datetick('x', 'dd-mmm-yyyy')
        xlabel('Date')
        ylabel('Difference from reference in %')
        title(sprintf('difference in percent for %s', sheet))
    else
        
        figure
        
        scatter(datetimes,difference_procent)
        xlabel('Date')
        ylabel('Difference from reference in %')
        title(sprintf('difference in percent for %s', sheet))
    end
end

%% Format dates
% This function takes in the data, the start date and the end date and then
% makes sure that the date is in the right format to be matched with the
% excell data. The function also looks for duplicates in the dates of the
% data since they need to be uniqe, if it finds a duplicate it returns them
% for the user to look at.


function dates_formated_unique = convertDates(data, start_date, end_date)
    % Lookup table
    monthNames = {'jan', 'feb', 'mars', 'apr', 'maj', 'juni', 'jul', 'aug', 'sep', 'okt', 'nov', 'dec','juli','mars'};
    englishMonthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec','Jul','Mar'};

    % Your dates
    dates = data(2,start_date:end_date);

    % Convert to datetime
    for i = 1:length(dates)
        dateString = lower(dates{i});  % Convert to lower case for case-insensitive comparison
        for j = 1:length(monthNames)
            dateString = strrep(dateString, monthNames{j}, lower(englishMonthNames{j}));  % Compare in lower case
        end
        dates{i} = datetime(dateString, 'InputFormat', 'dd-MMM-yyyy');  
    end
    dates_formated = [dates{:}];
    dates_formated_unique = unique(dates_formated);

    [C,~,idx] = unique(dates_formated,'stable'); % Get unique values and their indices
    n = accumarray(idx,1); % Count the occurrences
    duplicates = C(n>1); % Identify values with counts greater than 1

    if ~isempty(duplicates) || length(dates_formated) ~= length(dates_formated_unique)

        date_str = cellstr(datestr(duplicates)); 

        disp('There are mesurments done on the same date in the data, dates need to be unique. The duplicates are:')

        for i = 1 : length(date_str)

            sprintf('%s', date_str{i})

        end
    end
end

%% Number of columns
% This function finds the number of columns that containes data, this is
% becouse there are several cells that are surperflous and all we want is
% the data points.


function numNonNanInRow13 = number_of_columns(filename, sheetname)
    % Read Excel data into a table from a specific sheet
    data = readtable(filename, 'Sheet', sheetname);

    % Initialize counter
    numNonNanInRow13 = 0;

    % Loop over each column
    for i = 1:width(data)
        % Get cell in 13th row of current column
        cell = data{13, i};
        % If cell is numeric and not NaN, increment counter
        try
            if ~isnan(cell)
                numNonNanInRow13 = numNonNanInRow13 + 1;
            end
        catch
            % Ignore cells where isnan(cell) is not applicable
        end
    end
end

%% Select range
% This function asks the user for the data range that they want to analyse


function [start_date, end_date] = select_range()
    prompt = {'Enter start date (dd-mmm-yyyy):', 'Enter end date (dd-MMM-yyyy):'};
    dlgtitle = 'Input';
    date_range = inputdlg(prompt, dlgtitle);
    if isempty(date_range)
        disp('User selected Cancel')
    else
        start_date = date_range{1};
        end_date = date_range{2};
    end
end

%% Plot several energy levels
% This function plots all the data for the given energy levels in two plots
% one for 60 to 160 mev and one from 160 to 226 mev 

function diffrences_array = plot_several_energies(data,col_start, col_end, mev_several, datetimes, include_outliers, sheet)
    colorArray = [...
    0, 0, 0; % Black
    1, 0, 0; % Red
    0, 1, 0; % Green
    0, 0, 1; % Blue
    1, 1, 0; % Yellow
    0, 1, 1; % Cyan
    1, 0, 1; % Magenta
    0.5, 0.5, 0.5; % Gray
    0.8, 0.2, 0.2; % Dark Red
    0.2, 0.8, 0.2; % Dark Green
    0.2, 0.2, 0.8; % Dark Blue
    0.8, 0.8, 0.2; % Dark Yellow
    0.2, 0.8, 0.8; % Dark Cyan
    0.8, 0.2, 0.8; % Dark Magenta
    0.8, 0.4, 0.2; % Orange
    0.8, 0.6, 0.7; % Pink
    0.2, 0.4, 0.6; % Ocean Blue
    0.6, 0.4, 0.2; % Brown
    ];

    legendHandles1 = gobjects(0);
    legendNames1 = cell(0);

    legendHandles2 = gobjects(0);
    legendNames2 = cell(0);
    
    col_start = col_start - 1; 

    allowed_deviation = 10 ;
    
    diffrences_array = zeros(height(data) - 2, length(data) - 2);

    for mev_index = 1 : length(mev_several)

        mev =  mev_several(mev_index) ;
        index_mev = cellfun(@(x) isequal(x, mev), data(:,1));
        double_array = cell2mat(cellfun(@double, data(index_mev,2:end), 'UniformOutput', false));
        
        REF = double_array(1,1);

        if col_start == 2

            double_array = double_array(:, 2 : col_end - 1);
        else
            double_array = double_array(:, col_start : col_end);
        end
       
        % Calculate the number of columns in double_array
        n_cols = size(double_array, 2);  % -1 because we will discard the first column
        
        data_to_plot = zeros(1, n_cols);
        
        for i = 1 : n_cols 
            data_to_plot(1, i) = REF - double_array(1, i);
        end
        
        difference_procent = (data_to_plot ./ REF) * 100;
        diffrences_array(mev_index, : ) = difference_procent;


        if mev <= 160

            if include_outliers
                data_table = table(datetimes', difference_procent', 'VariableNames', {'Date', 'Data'});
                threshold = 10;
                is_above_threshold = data_table(data_table.Data > threshold, :) ;
                figure(10); % Switch to first figure
                h = scatter(is_above_threshold.Date, is_above_threshold.Data, 'filled', 'MarkerFaceColor', colorArray(mev_index, :));
                legendHandles1(end+1) = h;
                datetick('x', 'dd-mmm-yyyy')
                hold on;
            else
                difference_procent = replaceOutliers(difference_procent) ;
                figure(10); % Switch to first figure
                h = scatter(datetimes, difference_procent, 'filled', 'MarkerFaceColor', colorArray(mev_index, :));
                legendHandles1(end+1) = h;
                hold on;
                plot(datetimes, difference_procent, '--', 'Color', colorArray(mev_index, :));
                ylim([-1 7]);
                xlabel('Month')
                ylabel('Relative deviation from reference in %')
                title(sprintf('Relative deviation from reference in percent 60 to 160 MeV for %s', sheet))
                
            end
            legendNames1{end+1} = ['Energy ', num2str(mev), ' Mev'];

        elseif mev > 160

            if include_outliers
                data_table = table(datetimes', difference_procent', 'VariableNames', {'Date', 'Data'});
                threshold = 10;
                is_above_threshold = data_table(data_table.Data > threshold, :) ;
                figure(11); % Switch to second figure
                h = scatter(is_above_threshold.Date, is_above_threshold.Data, 'filled', 'MarkerFaceColor', colorArray(mev_index, :));
                legendHandles2(end+1) = h;
                datetick('x', 'dd-mmm-yyyy')
                hold on;
            else
                difference_procent = replaceOutliers(difference_procent) ;
                figure(11); % Switch to second figure
                h = scatter(datetimes, difference_procent, 'filled', 'MarkerFaceColor', colorArray(mev_index, :));
                legendHandles2(end+1) = h;
                hold on;
                plot(datetimes, difference_procent, '--', 'Color', colorArray(mev_index, :));
                ylim([-10 8]);
                xlabel('Month')
                ylabel('Relative deviation from reference in %')
                title(sprintf('Relative deviation from reference in percent 160 to 226 MeV for %s', sheet))
                
            end
            legendNames2{end+1} = ['Energy ', num2str(mev), ' Mev'];
        end
    end

    figure(10); % Switch to first figure
    line_1 = yline(allowed_deviation, '--r');
    line_1.LineWidth = 2;
    line_2 = yline(-allowed_deviation, '--r');
    line_2.LineWidth = 2;
    line_3 = yline(0, '--r');
    line_3.LineWidth = 1;

    if ~isempty(legendHandles1)
        legend(legendHandles1, legendNames1, 'Location', 'best');
    end
    hold off; 

    figure(11); % Switch to second figure
    line_4 = yline(allowed_deviation, '--r');
    line_4.LineWidth = 2;
    line_5 = yline(-allowed_deviation, '--r');
    line_5.LineWidth = 2;
    line_6 = yline(0, '--r');
    line_6.LineWidth = 1;

    if ~isempty(legendHandles2)
        legend(legendHandles2, legendNames2, 'Location', 'best');
    end
    hold off;
end

%% Takes out outliers 
% This function takes out the outliers from the data based on the outlier
% criterion set in the function


function data_no_outliers = replaceOutliers(data)

    outlierCriterion = 8 ;

    % Find the indices of the outliers
    outlierIndices = find(abs(data) > outlierCriterion);

    % Initialize a new array to hold the data with no outliers
    data_no_outliers = data;

    % For each outlier
    for i = 1:length(outlierIndices)

        % If the outlier is at the beginning of the array, replace it with the next data point
        if outlierIndices(i) == 1
            data_no_outliers(outlierIndices(i)) = data(2);

        % If the outlier is at the end of the array, replace it with the previous data point
        elseif outlierIndices(i) == length(data)
            data_no_outliers(outlierIndices(i)) = data(end-1);

        % If the outlier is in the middle of the array, replace it with the average of its two adjacent data points
        else
            data_no_outliers(outlierIndices(i)) = mean([data(outlierIndices(i)-1), data(outlierIndices(i)+1)]);
        end
    end
end

%% Linear regression and statistical analysis
% This function fits a line to the data of each energy level, the function
% then performes statistical analysis looking at the RMSE, R-value and
% P-value to determine if there is any statisticly significant change in
% the data.


function create_trend_lines(diffrences_array, sheet, mev_several)
    data = diffrences_array;
    n = height(data) - 1;
    p = zeros(n, 2);  % pre-allocating p for a 1st degree polynomial
    s = cell(n, 1);   % pre-allocating s
    errors = zeros(n, 1);  % pre-allocating errors
    R = zeros(n, 1);  % pre-allocating R values
    goodFit = false(n, 1);  % pre-allocating good fit check
    significant= false(n, 1);
    stdEner = zeros(1, length(n));
    pValue = zeros(1, length(n));

    % Energy levels array, modify this as per your data
    energyLevels = mev_several;  

    for energy = 1 : n
        ener = data(energy, 1:end);
        xx = 1:length(ener);
    
        [p(energy, :), s{energy}] = polyfit(xx, ener, 1);   
        [y, ~] = polyval(p(energy, :), xx, s{energy});  
        
        residuals = ener - y;  % residuals
        SSE = sum(residuals.^2);  % sum of squared residuals
        RMSE = sqrt(SSE / length(xx));  % root mean squared error
        
        errors(energy) = RMSE;
    
        % Calculate the correlation coefficient R
        correlationMatrix = corrcoef(xx, ener);
        R(energy) = correlationMatrix(1,2);  % The correlation coefficient is the off-diagonal element
        
        % Compare RMSE with the standard deviation of the energy level
        stdEner(energy) = std(ener);
        goodFit(energy) = RMSE < stdEner(energy);
    
        % Create a scatter plot for each energy level
        figure;
        plot(xx, ener, 'o', 'MarkerEdgeColor', 'blue'); % observed data
        hold on;
        plot(xx, y, 'LineWidth', 2, 'Color', 'red');  % Add the regression line
        hold off;
    
        title(sprintf('Energy Level %d: Observed vs Fitted for: %s', energyLevels(energy), sheet));
        xlabel('Time');
        ylabel('Relative deviation from reference in percent');
        legend('Observed', 'Fitted', 'Location', 'best');
        
        % p-value
        mdl = fitlm(xx, ener); % xx is your independent variable and ener is your dependent variable
        pValue(energy) = mdl.Coefficients.pValue(2); % pValue of the slope of the model
        
        significant(energy) = pValue(energy) < 0.05;

        % Add R, RMSE, std, slope, and p-value to the plot
        str = {sprintf('R: %.3f', R(energy)), sprintf('RMSE: %.3f', errors(energy)), sprintf('Std: %.3f', stdEner(energy)), sprintf('Slope: %.2f', p(energy, 1)), sprintf('p-value: %.5f', pValue(energy))};
        text(max(xx)*0.1, max(ener)*0.9, str, 'BackgroundColor', 'white');  % adjust position as necessary
        
    end
    
    % Plot RMSE and R as functions of energy level
    figure;
    plot(energyLevels, R, '-o', 'LineWidth', 1);
    title(sprintf('R as a function of Energy Level for: %s', sheet));
    xlabel('Energy Level (meV)');
    ylabel('R');
    yline(0, '--', 'Color', 'r')
    yline(1, '--', 'Color', 'r')
    yline(-1, '--', 'Color', 'r')
    hold off

    figure;
    plot1 = plot(energyLevels, errors, '-o','LineWidth', 1);  % errors is the RMSE
    hold on 
    plot2 = plot(energyLevels, stdEner, '-*', 'Color', 'r','LineWidth', 1);
    title(sprintf('RMSE and Std as a function of Energy Level for: %s', sheet));
    xlabel('Energy Level (meV)');
    ylabel('Value');
    legend([plot1 plot2], {'RMSE', 'Std'}, 'Location','northwest');
    hold off
    
    figure()
    plot(energyLevels, pValue, '-*', 'LineWidth', 1)
    hold on 
    yline(0.05 , '--', 'Color','r')
    xlabel('Energy level (meV)')
    ylabel('p-value')
    title(sprintf('P-value as a function of energy level for: %s', sheet))

end







