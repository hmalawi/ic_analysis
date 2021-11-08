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

    % Now call icrays3d.m to get the kernel-defining points
    % Need to store them in another cell array though, hope that works!
    pts{ii} = icrays3d(T, corelat, corelon, coredep, coredis, epid, ...
        mod, p, turnpt, 0)
end

end