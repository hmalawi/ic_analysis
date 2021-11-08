function iccoverage(fdir, fname, mod, vpahse, T)

% Open the file and read the data, skip the headerlines
% #Network, Station, sLatitude, sLongitude, EventID, tOrigin, eLatitude, eLongitude, Depth(km)
fid = fopen(strcat(fdir, fname), 'r');
data = textscan(fid, '%s%s%f%f%d%s%f%f%f', 'HeaderLine', 10);

% Now get the discretized ray for all of them




end