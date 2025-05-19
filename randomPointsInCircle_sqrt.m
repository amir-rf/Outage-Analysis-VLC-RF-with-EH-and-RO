function [x1, y1, x2, y2] = randomPointsInCircle_sqrt(R, numPoints)
    % R: radius of the circle
    % numPoints: number of random points to generate

    % Generate random angles for the first point
    theta1 = 2 * pi * rand(1, numPoints);

    % Generate random radii for the first point
    r1 = R * sqrt(rand(1, numPoints));

    % Convert polar coordinates to Cartesian coordinates for the first point
    x1 = r1 .* cos(theta1);
    y1 = r1 .* sin(theta1);

    % Generate random angles for the second point
    theta2 = 2 * pi * rand(1, numPoints);

    % Generate random radii for the second point
    r2 = R * sqrt(rand(1, numPoints));

    % Convert polar coordinates to Cartesian coordinates for the second point
    x2 = r2 .* cos(theta2);
    y2 = r2 .* sin(theta2);

    % Plot the circle and the generated points for visualization
%     figure;
%     hold on;
%     title('Random Points in a Circle');
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     axis equal;
%     viscircles([0, 0], R, 'EdgeColor', 'b', 'LineWidth', 2);
%     plot(x1, y1, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
%     plot(x2, y2, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'g');
%     hold off;
end