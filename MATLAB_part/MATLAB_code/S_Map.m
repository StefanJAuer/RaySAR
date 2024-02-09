function S_Map(Im_Ma,Bounce,a_p,r_p)

% function for displaying specular maps

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input: 
% - Im_Ma: matrix containing intensities
% - Bounce: string describing current bounce level
% - a_p: number of pixels in azimuth
% - r_p: number of pixels in range

% output: none

%--------------------------------------------------------------------------

% global parameters

global r_geom range_dir Output_path;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% Open new figure
figure;
set(gcf,'Name',Bounce,'Numbertitle','off')

% Range bottom up
if range_dir == 0
    Im_Ma = flipud(Im_Ma); % flip up down
    P = imagesc(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
else
    Im_Ma = fliplr(Im_Ma); % flip left right
    P = imagesc(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

colormap gray;

% display in ground range
if r_geom == 1
   xlabel('Azimuth Pixels'); ylabel('Ground Range Pixels'); title(Bounce);
end

% display in slant range
if r_geom == 0
   xlabel('Azimuth Pixels'); ylabel('Slant Range Pixels'); title(Bounce);
end

axis equal; axis image;

% Store image in folder
saveas(P,[Output_path,'/Maps/Figures/Specular/JPG/',Bounce,'.jpg'],'jpg');
saveas(P,[Output_path,'/Maps/Figures/Specular/FIG/',Bounce,'.fig'],'fig');
imwrite(Im_Ma,[Output_path,'/Maps/Frames/Specular/',Bounce,'.tif'],'tif','Compression','none');

