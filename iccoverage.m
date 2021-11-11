function iccoverage(fdir, fname, mod, vphase, T)

% Open the file and read the data, skip the headerlines
% #Network, Station, sLatitude, sLongitude, EventID, tOrigin, eLatitude, eLongitude, Depth(km)
fid = fopen(strcat(fdir, fname), 'r');
data = textscan(fid, '%s%s%f%f%d%s%f%f%f', 'HeaderLine', 10);

% Now get the discretized ray (with the kernels considered) for all of them
for ii = 1:length(data{1})
    % Call icray.m to get the descritized ray path
    [corelat, corelon, coredep, coredis, epid, p, turnpt, mod] = ... 
    icray(data{7}(ii), data{8}(ii), data{9}(ii), data{3}(ii), data{4}(ii),...
    mod, vphase);
   
    % Check whether the returned values are valid or not
    if isnan(corelat)
        continue
    else
        pts{ii} = icrays3d(T, corelat, corelon, coredep, coredis, epid, ...
            mod, p, turnpt, 0);
    end
end

% Remove empty cells from pts
pts = pts(~cellfun('isempty', pts));

% Now that we got all the points, choose some depths to check for
% First, find the ray that reaches the deepest (the one that has the most
% points.
[r, c] = cellfun(@size, pts);
% Store the depths somewhere as a reference (for easy access)
[val, index] = max(r);
dd = cell2mat(pts{index}(:,2));
% Maybe also store the maximum reached depth
[maxdep, maxpos] = max(dd);

ref = [dd(1), dd(21), dd(33), dd(46), dd(55), dd(71), dd(100), dd(122)];

% To store depths that already done
thisdone = zeros(length(ref),1);

% Choose some depths, or maybe do it for all? isn't that too much? maybe
% not! 'cause if didn't do it for all, how are you gonna choose?
for ii = 1:length(ref)
    % Check if this depth was already done
    findit = find(thisdone == ref(ii));
    if ~isempty(findit)
        continue
    end
    for jj = 1:length(pts)
        % Make temporary array of the depths of this patch for easy access
        temparr = cell2mat(pts{jj}(:,2));
        % Check if the current chosen depth exists
        pos = find(temparr == ref(ii));
        % If doesnt't exist or exists but empty, then skip this ray for now
        if ~exist('pos') || isempty(pos)
            continue
        end
        
        % Now get the sphere-cylinder intersection area for that depth
        % There will mostly be 2, so you need to check the length of pos
        for kk = 1:2
            % To avoid concatenation errors
            if length(pos) == 1 && kk == 2
                dim = size(v{1});
                v{jj,2} = zeros(dim(1), dim(2));
                break
            end
            % Assign the corresponding area points to an array
            pp = cell2mat(pts{jj}(pos(kk),1));
            % Need better arrangement to make things work
            pp = pp';
            
            % Call arrangepts.m to check teh points of the current polygon
            % (to avoid di
            % First to spherical coordinates
            [azi, ele, radi] = cart2sph(pp(:,1), pp(:,2), pp(:,3));
            % Then lat-lon
            lon = azi*180/pi;
            lat = ele*180/pi;
            % Now call inpolymoll to get proper matrices to plot
            [v{jj,kk}, xp, yp, xgr, ygr] = inpolymoll(lon, lat);
        end  
    end
    
    % Now need to add all the v we got
    S = sum(cat(3, v{:}), 3);
    % Plot them
    figure(1)
    clf
    pcolor(xgr, ygr, S); shading flat
    hold on
    plotcont([],[],2)
    plotplates([],[],2)
    hold off
    axis off image
    pstuff = sprintf('Maximum overlap %d', max(max(S)));
    text(1.65,-1.4,pstuff,'FontSize',6);
    % Save the figure maybe
    print(sprintf('/Users/hma/Desktop/maps/depth%d.png', ref(ii)),'-dpng','-r300')
    
    % Clear the variables S and v
    clear S
    clear v
    % Mark this depth as "done" to avoid repetition
    thisdone(ii) = ref(ii);
end

end