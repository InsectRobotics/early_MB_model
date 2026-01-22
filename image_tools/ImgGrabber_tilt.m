function img = ImgGrabber_tilt(x,y,z,th,X,Y,Z,colp,hfov,resolution)
% Grab an image from world model
% x,y,z : camera coordinates in metres
% th : heading direction in degrees, 0 is x+ axis
% X,Y,Z,colp : world data, grasses and color
% hfov : horizontal field of view in degrees, result image range from
%   [-hfov/2,hfov/2]
% resolution : image resolution in degrees/pixel
% Note: close figure when changing world parameters (X,Y,Z,colp), hfov and
%   resolution

z0=0;%getHeight(x,y,X,Y,Z);

%% prepare a figure or reuse figure in previous calls
f = 99;
if ~ishandle(f)
    figure(f);
    pause(0.25);
end
a = get(f,'CurrentAxes');
if isempty(a)
    figure(f);
    a = gca;
    set(f,'Color','c');
    set(f,'Name','Image Grabber');
%     set(f,'Renderer','OpenGL');
    set(f,'Renderer','zbuffer');
    imgwidth=ceil(hfov/resolution);
    imgheight = ceil(imgwidth/hfov*75);
    set(f, 'Position', [1 500 imgwidth imgheight]);
    set(a,'ActivePositionProperty','position');
    set(a, 'Position', [0 0 1 1]);
end

ph = findobj(get(a,'Children'),'UserData','plothere');
if isempty(ph)
    clear ph
else
    ph2 = findobj(get(a,'Children'),'UserData','plothere2');
    ph3 = findobj(get(a,'Children'),'UserData','plothere3');
    phg = findobj(get(a,'Children'),'UserData','plotground');    
end

%% transform world coordinates to camera coordinates
if(length(th) == 1) % fill tilt rotations if there is only a heading direction
    th = [th 0 0];
end
th = th/180*pi;
[dx,dy,dz] = tiltcamera(X-x,Y-y,abs(Z)-z-z0,th(1),th(2),th(3));

%% data projection from cartesian coordinates to spherical
[TH,PHI,R]=cart2sph(dx,dy,dz);

ind=(max(TH')-min(TH')<pi);
A1 = TH(ind,:);
E1 = PHI(ind,:);
D1 = R(ind,:);
c1 = colp(ind,:);

A2 = TH(~ind,:);
E2 = PHI(~ind,:);
D2 = R(~ind,:);
c2 = colp(~ind,:);

A3 = A2;
A3(A3<=0) = A3(A3<=0)+2*pi;
A4 = A2;
A4(A4>0) = A4(A4>0)-2*pi;

nn = 360;
Xp = 10*cos(linspace(-pi,pi,nn));
Yp = 10*sin(linspace(-pi,pi,nn));
%     Xp = [-10 -10 10.5 10.5]'; Yp = [-pi pi pi -pi]';
Zp = -z*ones(1,nn);
[Xp,Yp,Zp] = tiltcamera(Xp,Yp,Zp,th(1),th(2),th(3));
[Yp,Zp,Xp] = cart2sph(Xp,Yp,Zp);
[~,ip]=sort(Yp);
Xp = Xp(ip); Yp = Yp(ip); Zp = Zp(ip);
Xp = [0 Xp 0];
Yp = [Yp(1) Yp Yp(end)];
Zp = [-pi Zp -pi];
%% put data on the figure
if exist('ph','var')
    % change existing data
    set(ph,'XData',D1');
    set(ph,'YData',A1');
    set(ph,'ZData',E1');
    set(ph,'CData',c1');
    
    set(ph2,'XData',D2');
    set(ph2,'YData',A3');
    set(ph2,'ZData',E2');
    set(ph2,'CData',c2');
    
    set(ph3,'XData',D2');
    set(ph3,'YData',A4');
    set(ph3,'ZData',E2');
    set(ph3,'CData',c2');
    
    set(phg,'XData',Xp');
    set(phg,'YData',Yp');
    set(phg,'ZData',Zp');
    
    
else
    grasscolormap = zeros(64,3);
    grasscolormap(:,2) = linspace(0,1,64);
    % plot grasses
    ph = patch(D1',A1',E1',c1','EdgeColor','none');
    set(ph,'Userdata','plothere');
    colormap(grasscolormap);
    
    hold on
    ph2 = patch(D2',A3',E2',c2','EdgeColor','none');
    set(ph2,'Userdata','plothere2');
    colormap(grasscolormap);
    
    ph3 = patch(D2',A4',E2',c2','EdgeColor','none');
    set(ph3,'Userdata','plothere3');
    colormap(grasscolormap);

    groundcolor = [229 183 90]/255;
    phg = patch(Xp',Yp',Zp',groundcolor,'EdgeColor','none');
    set(phg,'Userdata','plotground');

    hold off
    
    % resize and crop
    % note: The view is looking in x+ axis direction.
    %  Thus, horizontal axis is y+ to the left
    %  or (hfov/2) to (-hfov/2) from left to right
    axis equal
    axis off
    hfov = hfov/180/2*pi;
    axis([0 14 -hfov hfov -pi/12 pi/3]);
    view([-90 0]);
end
drawnow
F = getframe(a);
img = F.cdata;


function x=pi2pi(x)
x=mod(x,2*pi);
x=x-(x>pi)*2*pi;

function [x1,y1,z1] = tiltcamera(x,y,z,yaw,pitch,roll)
s = size(x);
x = reshape(x,1,numel(x));
y = reshape(y,1,numel(y));
z = reshape(z,1,numel(z));
cz = cos(-yaw); sz = sin(-yaw);
cy = cos(-pitch); sy = sin(-pitch);
cx = cos(-roll); sx = sin(-roll);

v = ([1 0 0; 0 cx -sx; 0 sx cx]...
    *[cy 0 sy; 0 1 0; -sy 0 cy]...
    *[cz -sz 0; sz cz 0; 0 0 1]...
    )*[x;y;z];
x1 = reshape(v(1,:),s);
y1 = reshape(v(2,:),s);
z1 = reshape(v(3,:),s);
