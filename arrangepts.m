function [xnew, ynew, znew] = arrangepts(x, y, z)

% Find the distances
dist = sqrt([diff(x)].^2 + [diff(y)].^2 + [diff(z)].^2);
% Find the maximum distance, replace it with 0 in the original array to 
% avoid its contribution to the average. Then, find the average of the 
% rest of distances
[maxi, ind] = max(dist);
dist(ind) = 0;
avg = sum(dist)./(length(dist)-1);

% If that maximum did not differ too much from average (smooth transitions
% between the points)...
if maxi-avg < 0.1
    % Check the distance between the first and last point
    d = sqrt((x(1)-x(end)).^2 + (y(1)-y(end)).^2 + (z(1)-z(end)).^2);
    % If it's relatively large, fill the gap between them with some points
    % (linear interpolation). They should go from the last point to the
    % first point
    if d-avg >= 0.1
        % Define the points
        p1 = [x(end) y(end) z(end)];
        p2 = [x(1) y(1) z(1)];
        % Create, how about 100 points, between them? That should be enough
        n = 100;
        pset = linspace(0, 1, n)';
        pinterp = (1-pset)*p1 + pset*p2;
        % Exclude the first and last points 'cause they already exist in 
        % the original set of points
        pinterp = pinterp(2:end-1 , :);
        % Now retrun the new set of points (append the interpolated points)
        xnew = [x(:) ; pinterp(:,1)];
        ynew = [y(:) ; pinterp(:,2)];
        znew = [z(:) ; pinterp(:,3)];
    % If the distance is small, then there is no gap. Return the same
    % points
    elseif d-avg < 0.1
        xnew = x;
        ynew = y;
        znew = z;
    end
% If the maximum differs a lot from the average, then consider these points
% as the "gap" points
elseif maxi-avg >= 0.1
    % Define the points
    p1 = [x(ind) y(ind) z(ind)];
    p2 = [x(ind+1) y(ind+1) z(ind+1)];
    % Connect them
    n = 100;
    pset = linspace(0, 1, n)';
    pinterp = (1-pset)*p1 + pset*p2;
    % Return the new set of points
    xnew = [x(1:ind) ; pinterp(:,1) ; x(ind+1:end)];
    ynew = [y(1:ind) ; pinterp(:,2) ; y(ind+1:end)];
    znew = [z(1:ind) ; pinterp(:,3) ; z(ind+1:end)];
end

end