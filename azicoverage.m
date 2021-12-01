function azi = azicoverage(fdir, fname, pflag)
% azi = azicoverage(fdir, fname)
%
% This function is built to find and plot the azimuth of the rays
%
% INPUT:
%
% fdir           The directory at which the input file is located (this
%                should be in the same format of EQDATA output file)
% fname          Name of the data file
% pflag          1 to display plots
%
% OUTPUT:
%
% azi            An array holds all the azimuth values
%
%
% SEE ALSO:
%
% ICRAY
%
% Written by Huda Al Alawi (halawi@princeton.edu) - November 29, 2021
% Last modified by Huda Al Alawi - November 30, 2021
%

% Open the file and read the data, skip the headerlines
% #Network, Station, sLatitude, sLongitude, EventID, tOrigin, eLatitude, eLongitude, Depth(km)
fid = fopen(strcat(fdir, fname), 'r');
data = textscan(fid, '%s%s%f%f%d%s%f%f%f', 'HeaderLine', 10);

c = cell(length(data{1}), 4);

for ii = 1:length(data{1})
    % Call icray.m to get the descritized ray path in the inner core
    [corelat, corelon, coredep, coredis, epid, p, turnpt, mod] = ... 
    icray(data{7}(ii), data{8}(ii), data{9}(ii), data{3}(ii), data{4}(ii),...
    mod, vphase);

    % Store in- and out- inner core points
    % In-lat
    c{ii, 1} = corelat(1);
    % In-lon
    c{ii, 2} = corelon(1);
    % Out-lat
    c{ii, 3} = corelat(end);
    % Out-lon
    c{ii, 4} = corelon(end);
    
end

% Assign lat-lon coordinates of station and events to other variables for
% easy access
outlat = cell2mat(c(:,3));
outlon = cell2mat(c(:,4));
inlat = cell2mat(c(:,1));
inlon = cell2mat(c(:,2));

% Now define the x & y to then be able to convert into polar coordinates
y = sind(outlon-inlon) .* cosd(outlat);
x = (cosd(inlat) .* sind(outlat)) - (sind(inlat) .* cosd(outlat) .* cosd(outlon-inlon));

% Now, tan(theta) = y/x. Then, theta = atan(y/x)
theta = atand(y ./ x);

% Normalize the results to compass bearing - [0 360] range
azi = mod(theta+360, 360);

% Plot rose diagram
ax = polaraxes;
polarhistogram(azi, 'FaceColor', 'k')
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

end