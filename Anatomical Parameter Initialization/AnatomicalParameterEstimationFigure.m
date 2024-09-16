%% Anatomical Parameter Estimation Figure
clear all; close all; clc

mri_img = im2double(imread('..\Measurements_v1_no_lines.png'));
mri_img = mri_img(:,:,1)*(1/3) + mri_img(:,:,2)*(1/3) + mri_img(:,:,3)*(1/3);

[u,v] = size(mri_img);


points = readmatrix("AnatomicalMarkers.csv","NumHeaderLines",0);
points(:,2) = u-points(:,2);

I3_x = points(1:4,1);
I3_y = points(1:4,2);
I3 = polyshape(I3_x,I3_y);

DV_x = points(5:6,1);
DV_y = points(5:6,2);
AP_x = points(7:8,1);
AP_y = points(7:8,2);

R = mean([norm([diff(DV_x),diff(DV_y)]),norm([diff(AP_x),diff(AP_y)])])/2;



figure
imshow(flipud(mri_img),'XData',[0,v],'YData',[0,u]);
set(gca, 'YDir', 'normal');
set(gca,'FontName','Arial')
axis on;
hold all
plot(I3,"EdgeColor",0.5*[89,154,197]/255,"FaceColor",[89,154,197]/255,"FaceAlpha",0.3)
plot(DV_x,DV_y,"-r","LineWidth",2)
plot(AP_x,AP_y,"--r","LineWidth",2)



