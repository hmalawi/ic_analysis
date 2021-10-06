function pts = ...
         icrays3d(T, corelat, corelon, coredep, coredis, epid, mod, p, turnpt)
% pts = icrays3d(T, corelat, corelon, coredep, coredis, epid, mod, p, turnpt)
%
% This function returns a structure that has the kernel-defining-points
% of an inner core segment of a ray at different shells (depths).
%
% INPUT:
%
% T            The dominant period in seconds [Only for 1 and 2 seconds]
% corelat      The latitude of the discretized path of the innercore ray segment
% corelon      The longitude of the discretized path of the innercore ray segment
% coredep      The depth of the discretized path of the innercore ray segment
% raydis       The epicentral distance of the discretized path of the innercore ray segment [deg]
% epid         Epicentral distance
% mod          The chosen velocity model [defaulted]
% p            Ray parameter [s/deg]
% turnpt       Latitude, longitude, and depth of the turning point
%
% OUTPUT:
%
% pts         Cell array with all x,y,z points defining the "ray cylinders" in the first
%             column, the corresponding depth in the second column, and the
%             x,y,z points along the ray in the third column
%             these points are on a unit sphere of radius for the
%             corresponding depth 
%
% EXAMPLE:
%
% pts = icrays3d(T, corelat, corelon, coredep, coredis, epid, mod, p, ...
% turnpt);
% for ii = 1:length(pts)
%     pp = pts{ii,1};
%     c = pts{ii,3};
%     plot3(pp(1,:),pp(2,:),pp(3,:),'.-k')
%     hold on
%     plot3(c(1),c(2),c(3),'.-r', 'MarkerSize',10)
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     axis equal
%     axis([-1 1 -1 1 -1 1])
%     title(sprintf('%i',ii))
%     pause
%     hold on
% end
% hold off
%
% SEE ALSO:
%
% ICRAY
%
% Written by Huda Al Alawi - May 20th, 2021.
% Last modified by Huda Al Alawi - October 5th, 2021
%

% Define default values
defval('T', 1)
defval('mod', 'ak135')
defval('width', 300)

% We need to define the width of the kernels based on the dominant period
% and the epicentral distance (Calvet et al., 2006)
switch T
    case 1
        if epid>=146 && epid<160
            width = 295;
        elseif epid>=160 && epid<170
            width = 365;
        elseif epid>=170 && epid<180
            width = 390;
        end
    case 2
        if epid>=146 && epid<160
            width = 415;
        elseif epid>=160 && epid<170
            width = 510;
        elseif epid>=170 && epid<=180
            width = 540;
        end
 otherwise
  error('Need to specify an acceptable period')
end

% The radius of the inner core based on the chosen model...
a = eval(mod);
R = max(a.depth);
% The radius of the inner core will be R-(the first or last element of coredep array)
rsphere = R-coredep(1);

% A flag for the turning point
turn = 0;

% OF COURSE NEED TO INTERPOLATE THE DATA...

% The ray is discretized into small "cylinders". Each time, there will be
% in and out points. We should:
for ii = 1:length(corelat)-1
    
    % Before anything, check if we have reached the turning
    % point
    % If so, change the flag to 1, know it's index, define new arrays
    % for backward tracing, & skip this point for now...
    if turn==0 && corelat(ii)==turnpt(1) && corelon(ii)==turnpt(2) && ...
            coredep(ii)==turnpt(3)
        index = ii;
        % Change the flag so we don't mess things up
        turn = 1;
        % New arrays here
        newlat = flipud(corelat(index+1:end));
        newlon = flipud(corelon(index+1:end));
        newdep = flipud(coredep(index+1:end));
        % New counting variable
        jj = 1;
        continue
    end
    
    % Forward tracing until the turning point and then have to
    % trace it backward - it works!
    switch turn
        case 0
            % 1. Define these points, then call LINE3SPHERE to find the other
            % intersection point with the sphere of interest.
            % Be careful with the normalization process
            th = [corelon(ii), corelon(ii+1)] * pi/180;
            phi = [corelat(ii), corelat(ii+1)] * pi/180;
            r = (R-[coredep(ii), coredep(ii+1)])./rsphere;
            [x, y, z] = sph2cart(th(:), phi(:), r(:));
            xyz = line3sphere([x(1), y(1), z(1)], [x(2), y(2), z(2)], ...
                [0, 0, 0, (R-coredep(ii))/rsphere], 0);
            
            % Should find the coordinates of the second point to use CYLINDRIC
            % Make sure to select the right point
            if (x(1)==xyz(1,1)) && (y(1)==xyz(2,1)) && (z(1)==xyz(3,1))
                [th, phi, r] = cart2sph(xyz(1,2), xyz(2,2), xyz(3,2));
                outlon = th*180/pi; outlat = phi*180/pi;
            elseif x(1)==xyz(1,2) && y(1)==xyz(2,2) && z(1)==xyz(3,2)
                [th, phi, r] = cart2sph(xyz(1,1), xyz(2,1), xyz(3,1));
                outlon = th*180/pi; outlat = phi*180/pi;
            end
            
            % 2. Call CYLINDRIC to find the interection of a cylinder of radius r
            % (ray with kernel width) with a sphere. There will be top patch and bottom
            % patch. Take the top until the turning point, then will have to trace
            % it backwards
            % Radius of sphere and cylinder are normalized (using core's radius)
            [xyzS, topS, botS] = cylindric((width/2)/rsphere, [corelon(ii) corelat(ii)], ...
                [outlon outlat], (R-coredep(ii))/rsphere , 0);
            % Will take the upper patch
            prepts{ii,1} = topS;
            prepts{ii,2} = coredep(ii);
            % Store the x,y,z of the original point
            prepts{ii,3} = [x(1) y(1) z(1)];
            
        case 1
            % Back tracing using the new arrays, don't forget to replace
            % them ALL!
            % Find the other point using LINE3SPHERE
            th  =  [newlon(jj) newlon(jj+1)] * pi/180;
            phi =  [newlat(jj) newlat(jj+1)] * pi/180;
            r = (R-[newdep(jj) newdep(jj+1)])./rsphere;
            [x, y, z] = sph2cart(th(:), phi(:), r(:));
            xyz = line3sphere([x(1), y(1), z(1)], [x(2), y(2), z(2)], ...
                [0, 0, 0, (R-newdep(jj))/rsphere], 0);
            % Change coordinated and then use CYLINDRIC
            [th, phi, r] = cart2sph(xyz(1,2), xyz(2,2), xyz(3,2));
            outlon = th*180/pi; outlat = phi*180/pi;
            
            [xyzS, topS, botS] = cylindric((width/2)/rsphere, [newlon(jj) newlat(jj)], ...
                [outlon outlat], (R-newdep(jj))/rsphere, 0);
            % Will take the upper patch
            postpts{jj,1} = topS;
            postpts{jj,2} = newdep(jj);
            % Store the x,y,z of the original point
            postpts{jj,3} = [x(1) y(1) z(1)];
            
            % Update the counting variable jj
            jj = jj+1; 
    end
end

% Flip the points after turning point - just for better "animated" plotting
postpts = flipud(postpts);

% Now, combine prepts & postpts in one cell array
pts = [prepts; postpts];

end


