% Example 3D data
data = randn(100, 3);

% Perform PCA
coeff = pca(data);

% Create a 3D scatter plot of the data
scatter3(data(:, 1), data(:, 2), data(:, 3), '.');

hold on;

% Plot principal component vectors
origin = mean(data);  % Origin for quiver plot
scaleFactor = 1;  % Adjust this value to scale the length of the vectors

quiver3(0, 0,0, scaleFactor * coeff(1, 1), scaleFactor * coeff(2, 1), scaleFactor * coeff(3, 1), 'r', 'LineWidth', 2);
quiver3(0, 0, 0, scaleFactor * coeff(1, 2), scaleFactor * coeff(2, 2), scaleFactor * coeff(3, 2), 'g', 'LineWidth', 2);
quiver3(0,0, 0, scaleFactor * coeff(1, 3), scaleFactor * coeff(2, 3), scaleFactor * coeff(3, 3), 'b', 'LineWidth', 2);

hold off;

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Data with Principal Component Vectors');
grid on;