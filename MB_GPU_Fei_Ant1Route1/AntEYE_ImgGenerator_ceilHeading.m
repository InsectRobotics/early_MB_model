% generate data points and all the scanned images, stored in a struct.

load('./antview/AntData.mat');
load('./antview/world5000_gray.mat');

% retrieve a route from AntData
route = Ant1.InwardRouteData.Route1.One_cm_control_points/100;

% get an image every 10 cm
img_separation = 10; % [cm]
img_limit = floor(size(route,1)/10)*10; 
img_pos = route(1:img_separation:img_limit,:); % snapshot position
heading = zeros(size(img_pos,1)-1,1); % pointing angle of all the snapshot

% parameters for invoking ImgGrabber
eye_height = 0.01; % [m]
resolution = 4; % [degrees/pixel]
hfov = 296; % [degrees]
scan_range = 180; % [degrees]
scan_step = 2; % [degrees]


%% offset setup
num_pos = 80; % interval for 11 positions
num_start = 1; % starting position
%************************
% the idea is, when calculate the headings, we have to use two locations, 
% which means the first location in img_pos do not have a heading, so the 
% moving direction can only be applied to the second points. 
% Thus, I need 11 positions, to get 10 pos for training. 
% But I trained the net based only on the first 10 images before, which 
% actually was not precise enough. From now on just use the 11 pos.
%************************

num_off = 8; % 8 new positions corrsponding to each 10 pos.img_pos(1:10,:).
offDis = 0.02; % [m] distance of the offset positions
% off_pos = zeros(num_end, num_off); % 11 x 8 for now, only the last 10 has data

% save raw images in structure 
Raw_images = struct;
img_pos = img_pos(num_start:num_start+num_pos, :); % change here for different images

%% Get Headings for num_pos positions

for i = 1:num_pos
    heading(i) = atan2(img_pos(i+1,2)-img_pos(i,2),img_pos(i+1,1)-img_pos(i,1))*180/pi;
    % correspond to the second position.so start from 2
end

%% Round or ceil heading into recent even number:
for i = 1:num_pos
    if mod(floor(heading(i)), 2) == 0
        heading(i) = floor(heading(i));
    else
        heading(i) = ceil(heading(i));
    end
end

%% according to adjusted heading, update all the img_pos, except the first position.
step_size = 0.1; %[m]
for i = 1:num_pos
    img_pos(i+1,:) = img_pos(i,:) + [cos(heading(i)*pi/180), sin(heading(i)*pi/180)]*step_size;
end

%% Get all the corresponding images for all the positions in img_pos
% for i = 2:size(img_pos, 1) 
%     Raw_images.position(i).adjHeading = ImgGrabber(img_pos(i,1), img_pos(i,2), eye_height, heading(i), X,Y,Z,colp, hfov, resolution);
% end

%% Get all the scanned images for all the corresponding positions
% i.e., each position will have 180/2+1 = 91 images
for i = 1:num_pos
    for j = 1:scan_range/scan_step+1
        scan_angle = heading(i) + scan_range/2 - (j-1)*scan_step;
        Raw_images.position(i).scan_images(j).raw_image = ImgGrabber(img_pos(i,1), img_pos(i,2), eye_height, scan_angle, X, Y, Z, colp, hfov, resolution);
    end
end

%% Get the num_pos images as well as the offset ones
for i = 1:num_pos,
    for j = 1:num_off,
        temp_heading = heading(i)*pi/180 + (j-1)*(pi/4); % current heading in radians
        off_loc = img_pos(i,:) + offDis*[cos(temp_heading), sin(temp_heading)]; % cos-change in x, sin-change in y
        Raw_images.position(i).off_loc(j).raw_image = ImgGrabber(off_loc(1), off_loc(2), eye_height, heading(i), ... % all pointing to the same direc as the original pic
        X, Y, Z, colp, hfov, resolution);
        % add the original heading img into it
        Raw_images.position(i).origin = ImgGrabber(img_pos(i,1), img_pos(i,2), eye_height, heading(i), X, Y, Z, colp, hfov, resolution);
    end
end
%         Raw_images.position(i).off_pos(j).raw_image = ImgGrabber(img_pos(i,1)+,

%% Get 5 unknown testing images
num_test = 10;
route_test = Ant1.OutwardRouteData.Route1.One_cm_control_points/100;
img_pos_test = route_test(1:img_separation:img_limit,:);
img_pos_test = img_pos_test(num_start:num_start+num_test, :);
for i = 1:num_test,
    heading_test(i) = atan2(img_pos_test(i+1,2)-img_pos_test(i,2), img_pos_test(i+1,1)-img_pos_test(i,1))*180/pi;
    Raw_images.test_position(i).origin = ImgGrabber(img_pos_test(i,1), img_pos_test(i,2), eye_height, heading_test(i), X, Y, Z, colp, hfov, resolution);
end
Raw_images.img_pos = img_pos;
save(sprintf('Raw_images_numPos_180degree%d_win_evenHeading.mat',num_pos), 'Raw_images');
        

    