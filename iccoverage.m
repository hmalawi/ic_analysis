function S = iccoverage(fdir, fname, mod, vphase, T, outdir, dmin, dmax, pflag)
% S = iccoverage(fdir, fname, mod, vphase, T, outdir, dmax, pflag)
%
% This function is built to plot inner-core spatial coverage maps.
%
% INPUT:
%
% fdir           The directory at which the input file is located (this
%                should be in the same format of EQDATA output file)
% fname          Name of the data file
% mod            The chosen velocity model [defaulted]
% vphase         Seismic velocity phase (i.e., PKIKP for the innercore) [defaulted]
% T              The dominant period [defaulted]	
% outdir         The directory at which the maps will be saved
% dmin           The index of starting data point [defaulted: 1]
% dmax           How far in the data file to progress [default: all]
% pflag          1 to save plots
%
% OUTPUT:
%
% S              A cell array that has plottable Mollweide heat maps at
%                different depths of the inner-core in the first column,
%                the corresponding depth in the secon column, and the
%                maximum overlap in the third column
%
% In addition to a set of maps will be generated and saved into "outdir"
%
%
% SEE ALSO:
%
% ICRAY, ICRAYS3D
%
% Written by Huda Al Alawi (halawi@princeton.edu) - November 6, 2021
% Last modified by Huda Al Alawi - November 30, 2021
%

% Open the file and read the data, skip the headerlines
% #Network, Station, sLatitude, sLongitude, EventID, tOrigin, eLatitude, eLongitude, Depth(km)
fid = fopen(strcat(fdir, fname), 'r');
data = textscan(fid, '%s%s%f%f%d%s%f%f%f', 'HeaderLine', 10);

defval('dmax', length(data{1}))
defval('dmin', 1)
defval('pflag', 0)

% Initialize PTS with the max possible number of points later

% Independent counting index for the 
jj = 1;

% Now get the discretized ray (with the kernels considered) for all of them
for ii = dmin:dmax
    % Call icray.m to get the descritized ray path
    [corelat, corelon, coredep, coredis, epid, p, turnpt, mod] = ... 
    icray(data{7}(ii), data{8}(ii), data{9}(ii), data{3}(ii), data{4}(ii),...
    mod, vphase);
   
    % Check whether the returned values are valid or not
    if isnan(corelat)
        % Update the index and continue
        jj = jj + 1;
        continue
    else
        % Save and update the index
        pts{jj} = icrays3d(T, corelat, corelon, coredep, coredis, epid, ...
            mod, p, turnpt, 0);
        jj = jj + 1;
    end
    % Percentage update
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

% Choose some depths, you can do it all though
refs=[1 21 33 46 55 71 100 maxpos];
ref =dd(refs(refs<=length(dd)));

% To store depths that are already done
thisdone = zeros(length(ref),1);

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
            
            % Call arrangepts.m to check the points of the current polygon
            % (to avoid distortions when plotting)
            [x, y, z] = arrangepts(pp(:,1), pp(:,2), pp(:,3));
            
            % First to spherical coordinates
            [azi, ele, radi] = cart2sph(x, y, z);
            % Then lat-lon
            lon = azi*180/pi;
            lat = ele*180/pi;
            % Now call inpolymoll to get proper matrices to plot. Store
            % them all into a cell array
            [v{jj,kk}, xp, yp, xgr, ygr] = inpolymoll(lon, lat);
        end  
    end
    
    % Now need to add all the v(s) we got for this depth and store it 
    % in another cell array (I could have generated the maps here.
    % However, since the colorbar range should be uniform for all depths,
    % I need to know the maximum overlap at all depths before deciding.)
    S{ii,1} = sum(cat(3, v{:}), 3);
    S{ii,2} = ref(ii);
    S{ii,3} = max(max(S{ii,1}));
    
    % Mark this depth as "done" to avoid repetition
    thisdone(ii) = ref(ii);
    % Clear the variable S v to make sure it won't contribute to the
    % other depths
    clear v

end

% When done with summing the points, save "S" in outdir. Remember to
% indicate the chosen indices in the name
save(sprintf('sum%d-%d', dmin, dmax), 'S');


if pflag == 1
    % Find the maximum and minimum overlap for uniform colorbar range
    maxx = max([S{:,3}]);
    minn = min([S{:,3}]);
    
    % Save the heat maps without displaying
    for h = 1:length(ref)
        outname = sprintf('%sdepth%d.png', outdir, ref(h));
        f = figure('visible', 'off');
        clf
        pcolor(xgr, ygr, S{h,1}); shading flat
        hold on
        plotcont([], [], 2)
        plotplates([], [], 2)
        colorbar
        %caxis([min(0,minn) maxx])
        caxis([minn maxx])
        hold off
        axis off image
        pstuff = sprintf('The maximum overlap %d', S{h,3});
        text(1.65, -1.4, pstuff, 'FontSize', 6);
        print(outname, '-dpng', '-r300')
    end
end

end