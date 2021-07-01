function [corelat,corelon,coredep,epid,p,turnpt] = ...
    icray(eqlat,eqlon,eqdepth,stalat,stalon,mod,vphase)
% [corelat,corelon,coredep,epid] = ...
%                      icray(eqlat,eqlon,eqdepth,stalat,stalon,mod,vphase);
%
% This function is mainly built to calculate find the the discretized path
% of the inner core segment of a ray.
%
% INPUT:
%
% eqlat           Earthquake's latitude
% eqlon           Earthquake's longitude
% stalat          Station's latitude
% stalon          Station's longitude	
% mod             The chosen velocity model [defaulted]
% vphase          Seismic phase (i.e., PKIKP for the innercore) [defaulted]
%
% OUTPUT:
%
% corelat         The latitude of the discretized path of the innercore ray segment
% corelon         The longitude of the discretized path of the innercore ray segment
% coredep         The depth of the discretized path of the innercore ray segment
% turnpt          The turning point (coordinates and depth)
% epid            Epicentral distance
% p               Ray parameter [s/deg]
%
%
% Written by Huda Al Alawi - May 16th, 2021.

% Define default values
defval('mod','ak135')
defval('vphase','PKIKP')

% Now let's use taup to discretize some parameters
data=tauppath('mod',mod,'dep',eqdepth,'ph',vphase,'sta',[stalat stalon],...
     'evt',[eqlat eqlon]);
 
% Epicentral distance
epid=data.distance;
% Ray parameter
p=data.rayparameter;
% Discretize epicentral distance
raydis=data.path.distance;
% Discretize depth
raydep=data.path.depth;
% Discretize location
raylat=data.path.latitude;
raylon=data.path.longitude;

% Let's find the parameters of the selected Earth model
a=eval(mod);
% What's the radius of the Earth?
R=max(a.depth);
% To find the angle of incidence using rayparameter, we need to know the
% velocity at each point, so...
for index=1:length(raydep)
  eval(sprintf('vp(index)=%s',sprintf('getfield(%s(''depths'',%f),''vp'');'...
      ,mod,raydep(index))));
end
% Incidence angle - protected against rounding imaginary numbers
% check REALIZE for other cases not relevant here
iang=asin(min(1,p*vp(:)./[R-raydep(:)]*180/pi))*180/pi;
% This is because, for a value that is ~ 1 (at 90 degree, mostly turning 
% point), will start to get complex numbers.

% return the turning point info.
iturn=find(iang==90);
turnpt=[raylat(iturn),raylon(iturn),raydep(iturn)];

% Now, we only need the discretized path for the inner core segment
% What is the depth of the inner core?
index=find(a.vs==0);
icdep=a.depth(index(end)+1);

% Data of innercore only...
inout=find(raydep==icdep);
% If, for any reason, there is a problem in the ray, i.e., not having two
% points at the depth of the inner core, return nothing
if length(inout)~=2
    corelat=NaN; corelon=NaN; coredep=NaN; turnpt=NaN; epid=NaN; p=NaN;
    return
% If everything went well, return innercore data
else
    corelat=[raylat(inout(1):inout(2))];
    corelon=[raylon(inout(1):inout(2))];
    coredep=[raydep(inout(1):inout(2))];
    
end


end
