function arrangepts(x, y, z)

% Find the distances
dist = sqrt([diff(x)].^2 + [diff(y)].^2 + [diff(z)].^2);
% Find the maximum distance, replace it with 0 in the original array, then
% find the average of the rest of distances
[maxi, ind] = max(dist);
dist(ind) = 0;
avg = sum(dist)./(length(dist)-1);

% If that maximum did not differ too much from average (smooth transitions
% between the points)...
if maxi-avg < 0.1
    % Check the distance between the first and last point
    d = sqrt((x(1)-x(end)).^2 + (y(1)-y(end)).^2 + (z(1)-z(end)).^2);
    % If it's relatively large, fill the gap between them with some points
    % (linear interpolation)
    
    
    
    
    
end

end