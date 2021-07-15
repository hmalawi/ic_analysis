function [pts,index]=icrays3d(T,corelat,corelon,coredep,epid,mod,p,turnpt)
% pts=icrays3d(T,corelat,corelon,coredep,epid,mod,p,turnpt);
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
% epid         Epicentral distance
% mod          The chosen velocity model [defaulted]
% p            Ray parameter [s/deg]
% turnpt       Latitude,longitude, and depth of the turning point
%
% SEE ALSO:
%
% ICRAY
%
% Written by Huda Al Alawi - May 20th, 2021.
% Last modified by Huda Al Alawi - July 2nd, 2021
%

% Define default values
defval('mod','ak135')

% We need to define the width of the kernels based on the dominant period
% and the epicentral distance (Calvet et al., 2006)
switch T
    case 1
        if epid>=146 && epid<160
            width=295;
        elseif epid>=160 && epid<170
            width=365;
        elseif epid>=170 && epid<180
            width=390;
        end
       
    case 2
        if epid>=146 && epid<160
            width=415;
        elseif epid>=160 && epid<170
            width=510;
        elseif epid>=170 && epid<180
            width=540;
        end
end

% The radius of the inner core based on the choosen model...
a=eval(mod);
R=max(a.depth);
% The radius of the inner core will be R-(the first or last element of coredep array)
rsphere=R-coredep(1);

% A flag for the turning point
turn=0;

% OF COURSE NEED TO INTERPOLATE THE DATA...

% The ray is discretized into small "cylinders". Each time, there will be
% in and out points. We should:
for ii=1:length(corelat)-1
    % 1. Define these points, then call LINE3SPHERE to find the other
    % intersection point with the sphere of interest (specific depth of the inner core).
    th=[corelon(ii), corelon(ii+1)]*pi/180;
    phi=[corelat(ii), corelat(ii+1)]*pi/180;
    r=(R-[coredep(ii), coredep(ii+1)])./R;
    [x,y,z]=sph2cart(th(:),phi(:),r(:));
    xyz=line3sphere([x(1),y(1),z(1)],[x(2),y(2),z(2)],[0,0,0,(R-coredep(ii))/R],0);
    
    % Should find the coordinates of the second point to use CYLINDRIC
    [th,phi,r]=cart2sph(xyz(1,2),xyz(2,2),xyz(3,2));
    outlon=th*180/pi; outlat=phi*180/pi;
    
    % 2. Call CYLINDRIC to find the interection of a cylinder of radius r
    % (ray with kernel width) with a sphere. There will be top patch and bottom
    % patch. Take the top until the turning point, then take the bottom.
    % Radius of sphere and cylinder are normalized (using sphere's radius)
    [xyzS,topS,botS]=cylindric((width/2)/rsphere,[corelon(ii) corelat(ii)],[outlon outlat],1,0);
    
    % Check if have reached the turning point yet
    % If so, change the flag to 1;
    if corelat(ii)==turnpt(1) && corelon(ii)==turnpt(2) && coredep(ii)==turnpt(3)
        turn=1;
        index=ii;
    end
    
    % Top or bottom?
    % Check, turn the points and the depths
    if turn==0
        pts{ii,1}=topS;
        pts{ii,2}=coredep(ii);
        
    elseif turn==1
        pts{ii,1}=botS;
        pts{ii,2}=coredep(ii+1);
    end

end


end
