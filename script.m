clear all;close all;clc;
format long;

%% kml read
% Different functions to get list of coordinates from kml polygon
% Polygon is generated using the drawing tools on https://map.geo.admin.ch
% 
% Attention: The exportet *.kml ftom geo.admin is not formated correctly
% (everything on one line), therefore all kml importing functions in MATLAB
% don't work. Workaround: Open .kml in Googel Earth and export it again.



% map_structure = kml_shapefile('Placemark.kml');
% [lat,lon,z] = read_kml('Placemark.kml');

% map_structure = kml2struct('Placemark.kml'); %kml of zurichsee

map_structure = kml2struct('Technosee.kml'); %kml of technosee


%% interpolate shorline points to no bigger distance than X meters
maxdiff_m= 10;
maxdiff=distdim(maxdiff_m,'meters','degrees');

[latout,lonout] = interpm(map_structure.Lat,map_structure.Lon,maxdiff);
latout=latout(1:end-1);
lonout=lonout(1:end-1);

%% calc shoreline 150m
%bufferm returns not only the new line but a polygon which contains the
%whole bufferzone.
%Do not calculate bufferm from interpolated polynom, takes too long

dist_m= 150;
dist_deg=distdim(dist_m,'meters','degrees');

[latb,lonb] = bufferm(map_structure.Lat,map_structure.Lon,dist_deg,'in');

%% crop inner polygon
%the direction of the polygon is not every time the same. Check plot if the
%cropping is done right

%find where polygon is connected
index = find(latb==latb(end),1,'first');
index2 = find(lonb==lonb(end),1,'first');

% cut #1
% latb=latb(index:end);
% lonb=lonb(index2:end);

% cut #2
latb=latb(1:index);
lonb=lonb(1:index2);

% remove NaN
latb=latb(1:end-2);
lonb=lonb(1:end-2);



%% calc shoreline 300m
dist_m2= 300;
dist_deg2=distdim(dist_m2,'meters','degrees');

[latb2,lonb2] = bufferm(map_structure.Lat,map_structure.Lon,dist_deg2,'in');


%% crop inner polygon

index3 = find(latb2==latb2(end),1,'first');
index4 = find(lonb2==lonb2(end),1,'first');

% latb2=latb2(index3:end);
% lonb2=lonb2(index4:end);
latb2=latb2(1:index3);
lonb2=lonb2(1:index4);

latb2=latb2(1:end-2);
lonb2=lonb2(1:end-2);

%% plot maps
% should be 3 concentric polygons

% geoshow(map_structure)
figure(10);
hold on
geoshow(latb,lonb)
% geoshow(latbi,lonbi)
geoshow(latb2,lonb2)
geoshow(latout,lonout)

grid on


%% write files

% interpolated and exported below: % dlmwrite('zurichsee_border1.csv', [latb,lonb], 'delimiter', ',', 'precision', 9); 
% dlmwrite('zurichsee_border2.csv', [latb2,lonb2], 'delimiter', ',', 'precision', 9); 
% dlmwrite('zurichsee_shore.csv', [latout,lonout], 'delimiter', ',', 'precision', 9); 

% interpolated and exported below: dlmwrite('technosee_border1.csv', [latb,lonb], 'delimiter', ',', 'precision', 9); 
dlmwrite('technosee_border2.csv', [latb2,lonb2], 'delimiter', ',', 'precision', 9); 
dlmwrite('technosee_shore.csv', [latout,lonout], 'delimiter', ',', 'precision', 9); 

% csvwrite('csvlist.csv',[latb,lonb])%nur 4 kommastellen präzision



% kmlwriteline('filename.kml',latb,lonb);
% kmlwriteline('filename2.kml',latb2,lonb2);

if 0
%% latlong to swiss

[xb1,yb1] = deg2ch1903plus(latout,lonout);
[xb2,yb2] = deg2ch1903plus(latb,lonb);
[xb3,yb3] = deg2ch1903plus(latb2,lonb2);


%%

figure(20)
hold on
plot(xb1,yb1,'b')
plot(xb2,yb2,'r')
plot(xb3,yb3,'g')


%% grid

shiftx= 682000;
xb1=xb1-shiftx;
xb2=xb2-shiftx;
xb3=xb3-shiftx;

shifty= 228000;
yb1=yb1-shifty;
yb2=yb2-shifty;
yb3=yb3-shifty;


figure(21);hold on;


plot(xb1,yb1,'b')
plot(xb2,yb2,'r')
plot(xb3,yb3,'g')


%%

% grid2 = 1000*ones(20000/10,25000/10);


[iv,jv]=ndgrid(1:10:25000,1:10:20000);


INb1 = inpolygon(iv,jv,xb1,yb1);
INb2 = inpolygon(iv,jv,xb2,yb2);
INb3 = inpolygon(iv,jv,xb3,yb3);


in_matrix=1*ones(size(iv));
in_matrix(INb1)=1/50;
in_matrix(INb2)=1/20;
in_matrix(INb3)=1/10;

surf(in_matrix)
end %code no longer used

%% 

%latlong border
latb1=latout;
lonb1=lonout;

%latlong 300
latb3=latb2;
lonb3=lonb2;
maxdiff_m= 10;
maxdiff=distdim(maxdiff_m,'meters','degrees');

[latb3,lonb3] = interpm(latb3,lonb3,maxdiff);

%latlong 150
latb2=latb;
lonb2=lonb;
[latb2,lonb2] = interpm(latb2,lonb2,maxdiff);


% dlmwrite('zurichsee_border1.csv', [latb2,lonb2], 'delimiter', ',', 'precision', 9); 
dlmwrite('technosee_border1.csv', [latb2,lonb2], 'delimiter', ',', 'precision', 9); 

% dlmwrite('zurichsee_border2.csv', [latb2,lonb2], 'delimiter', ',', 'precision', 9); 
% dlmwrite('zurichsee_shore.csv', [latout,lonout], 'delimiter', ',', 'precision', 9); 


%% 
figure(30);
hold on
geoshow(latb1,lonb1)
geoshow(latb2,lonb2)
geoshow(latb3,lonb3)

%point on shore
pshore_1_lat=47.318850918570940;
pshore_1_lon=8.553807722787102;
geoshow(pshore_1_lat,pshore_1_lon,'DisplayType','point')

%point on 150
p150_1_lat=47.319106029904590;
p150_1_lon=8.555897406912782;
geoshow(p150_1_lat,p150_1_lon,'DisplayType','point')




%%




for i=1:length(latb1)
[arclen1(i),az] = distance(pshore_1_lat,pshore_1_lon,latb1(i),lonb1(i));
dist_meter1(i)=distdim(arclen1(i),'degrees','meters');
[arclen2(i),az] = distance(p150_1_lat,p150_1_lon,latb1(i),lonb1(i));
dist_meter2(i)=distdim(arclen2(i),'degrees','meters');
if mod(i,100)
    i
end
end


figure
plot(dist_meter1)
hold on
plot(dist_meter2)













