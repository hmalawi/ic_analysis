function S = azicoverage(fdir, fname)

% Open the file and read the data, skip the headerlines
% #Network, Station, sLatitude, sLongitude, EventID, tOrigin, eLatitude, eLongitude, Depth(km)
fid = fopen(strcat(fdir, fname), 'r');
data = textscan(fid, '%s%s%f%f%d%s%f%f%f', 'HeaderLine', 10);

% Assign lat-lon coordinates of station and events to other variables for
% easy access
slat = data{3};
slon = data{4};
elat = data{7};
elon = data{8};

% Now define the x & y to then be able to convert into polar coordinates
y = sind(slon-elon) .* cosd(slat);
x = (cosd(elat) .* sind(slat)) - (sind(elat) .* cosd(slat) .* cosd(slon-elon));

% Now, tan(theta) = y/x. Then, theta = atan(y/x)
theta = atand(y ./ x);

% Normalize the results to compass bearing - [0 360] range
azi = mod(theta+360, 360);



end