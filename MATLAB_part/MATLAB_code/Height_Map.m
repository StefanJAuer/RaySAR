function Height_Map(Im_Ma,a_p,r_p,H_Text)

% function for displaying height map (maximum elevation coordinate)

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input: 
% - Im_Ma: matrix containing maximum elevation heights (for each pixel)
% - a_p: number of pixels in azimuth
% - r_p: number of pixels in range
% - H_Text: plot title (string)

% output:
% - Figure stored in Maps/Figures/Elevation\_Map/, image stored in
% Maps/Frames/Elevation\_Map/
% - values of image matrix stored in Maps/Figures/Elevation\_Map/

%--------------------------------------------------------------------------

% global parameters

global r_geom range_dir ang Output_path;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% 1.) Save values of image matrix to text file

% remove old version
if exist([Output_path,'/Maps/Figures/Elevation_Map/El_Map.txt'],'file') == 2
   
    % delete
    delete([Output_path,'/Maps/Figures/Elevation_Map/El_Map.txt']);
end

% save new version
save([Output_path,'/Maps/Figures/Elevation_Map/El_Map.txt'],'Im_Ma','-ASCII');

%--------------------------------------------------------------------------

% 2.) Display elevation map

min_Im = min(min(Im_Ma));
max_Im = max(max(Im_Ma));
  
% Scale image entries to interval of gray values
num_g = 2^8; % maximum gray value
Im_Ma = im_trans(Im_Ma,min_Im,max_Im,num_g-1);

% Open new figure
figure;
set(gcf,'Name',H_Text,'Numbertitle','off')

if range_dir == 0
    Im_Ma = flipud(Im_Ma); % flip up down
    P = image(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

if range_dir == 1
    Im_Ma = fliplr(Im_Ma); % flip left right
    P = image(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

% Colormap --> gray values
Map = colormap(jet(num_g));

% Colorbar
colorbar;

% display in ground range 
if r_geom == 1
    xlabel('Azimuth Pixels'); ylabel('Ground Range Pixels'); title([H_Text ' (dB) - angle of incidence: ' num2str(ang*(180/pi)) ' degrees']);
end

% display in slant range
if r_geom == 0
    xlabel('Azimuth Pixels'); ylabel('Slant Range Pixels'); title([H_Text ' (dB)']);
end

axis equal; axis image;

% Store image in folder
saveas(P,[Output_path,'/Maps/Figures/Elevation_Map/',H_Text,'.jpg'],'jpg');
saveas(P,[Output_path,'/Maps/Figures/Elevation_Map/',H_Text,'.fig'],'fig');

%--------------------------------------------------------------------------