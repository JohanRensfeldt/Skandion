% 1. Define the Gaussian function
gaussian = @(x, mu, sigma) 0.7 * exp(-0.5 * ((x - mu) / sigma).^2); % 0.7 amplitude to make the peak 70%

% 2. Generate x-values
x = -10:0.01:10; % Span to accommodate 15 gaussians.

% 3. Define means, standard deviation
numGaussians = 20;
spacing = 0.5; % Spacing between means
startPoint = -spacing * (numGaussians-1)/2; % To center the Gaussians around 0.
mus = startPoint:spacing:(startPoint + spacing*(numGaussians-1));
original_mus = mus; % Saving original positions

% Apply the displacement for two of the Gaussians. 
mus(7) = mus(7) - 0.15; % -1.5 mm 
%mus(9) = mus(9) + 0.15; % +1.5 mm

sigma = 0.3; % 3 mm

% Compute the Gaussians and their sum for the displaced Gaussians
allGaussians_displaced = zeros(length(mus), length(x));
for i = 1:length(mus)
    allGaussians_displaced(i, :) = gaussian(x, mus(i), sigma);
end
sumGaussians_displaced = sum(allGaussians_displaced, 1);

% Compute the Gaussians and their sum for the original Gaussians (without displacement)
allGaussians_original = zeros(length(original_mus), length(x));
for i = 1:length(original_mus)
    allGaussians_original(i, :) = gaussian(x, original_mus(i), sigma);
end
sumGaussians_original = sum(allGaussians_original, 1);

% Normalize both sums to the peak of the original (non-displaced) sum
normalization_factor = max(sumGaussians_original);
sumGaussians_displaced = sumGaussians_displaced / normalization_factor;
sumGaussians_original = sumGaussians_original / normalization_factor;


% Plot
figure;
hold on;

% Plot the individual displaced Gaussians
for i = 1:length(mus)
    plot(x, allGaussians_displaced(i, :), 'LineStyle', '--', 'Color', [0.8 0.8 0.8]); % Plot in light gray
end

% Plot the normalized sum of displaced Gaussians
plot(x, sumGaussians_displaced, 'LineWidth', 2, 'Color', 'red');

% Plot the normalized sum of original Gaussians (without displacement)
plot(x, sumGaussians_original, 'LineWidth', 2, 'Color', 'blue');

% Calculating dose error
doseError = max(sumGaussians_displaced) - max(sumGaussians_original);

% Displaying the information on the plot
infoString = sprintf('Dose Error: %0.2f%%\nSpot Size (sigma): %0.2f cm\nX Displacement: +1.5 mm', doseError*100, sigma);
textPos = [min(x) + 0.2, 0.9];
text(textPos(1), textPos(2), infoString, 'FontSize', 10, 'BackgroundColor', 'white');

xlabel('Position (cm)');
ylabel('Relative Dose %');
title('Comparison of Sum of Gaussians with and without Displacements');
hold off;