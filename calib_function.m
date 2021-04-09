function [T]=calib_function(calibpoints)

%%%
%% 1 :  350 412
%% 2 :  350 662
%% 3 :  400 662	s
%%
if ~nargin
    calibpoints = [412 350 1; 662 350 1 ; 662 400 1]';

else
    calibpoints = calibpoints;
end

supportedtypes={'*.png','PNG Files'; '*.bmp','BMP Files'};

currentFolder = pwd;
[fname,pname,typeind]=uigetfile(supportedtypes,'Choose Calibration image',currentFolder);

info = imfinfo([pname,fname]);

height= info(1).Height;
width=info(1).Width;

F = imread(fullfile(pname, fname));

[m,n,o] = size(F);      
figure;

img = imagesc((double(F)));
colormap(gray);
title('Choose reference points in right order');
pts=ginput_modified(3);

x=round(pts(:,1)); y=round(pts(:,2));
pause(0.7);
d=30;
for k=1:3
    set(gca,'ylim',[y(k)-d y(k)+d]);
    set(gca,'xlim',[x(k)-d x(k)+d]);
    title(sprintf('Point #%d - zoom',k));
    pts1=ginput_modified(1);
    x1(k)=pts1(1,1); 
    y1(k)=pts1(1,2);
end
campoints = [x1; y1; [1, 1, 1]];                                                                                                                                                                                                                                                         

T = calibpoints*inv(campoints);

T(3,:) = [0 0 1];
sintheta = 0.5*(-T(1,2) + T(2,1));
T (1,2) = -sintheta;
T(2,1) = sintheta;

close;


