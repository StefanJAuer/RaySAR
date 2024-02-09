function Histo_Map(Im_Ma,a_p,r_p,scat_max,H_Text)

% function for displaying 2D-histogram 
% --> number of scatterers per resolution cell

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input: 
% - Im_Ma: matrix containing number of scatterers for each resolution cell
% - a_p: number of pixels in azimuth
% - r_p: number of pixels in range
% - scat_max: maximum number of scatterers within matrix Im_Ma (used for displaying purposes)
% - H_Text: string for the title of the plot

% output:
% - image stored in Maps/Figures/Scatterer_Map/

%--------------------------------------------------------------------------

% global parameters

global r_geom range_dir Output_path;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% Open new figure
figure;
set(gcf,'Name',H_Text,'Numbertitle','off')

if range_dir == 0
    Im_Ma = flipud(Im_Ma); % flip up down
    P = imagesc(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

if range_dir == 1
    Im_Ma = fliplr(Im_Ma); % flip left right
    P = imagesc(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

% display in ground range 
if r_geom == 1
    xlabel('Azimuth Pixels'); ylabel('Ground Range Pixels'); title(H_Text);
end

% display in slant range
if r_geom == 0
    xlabel('Azimuth Pixels'); ylabel('Slant Range Pixels'); title(H_Text);
end

axis equal; axis image;

% Color settings
Map = colormap(jet(scat_max+1));
color_set = colorbar;
set(color_set,'YTick',0:1:scat_max);
set(color_set,'TickLength',0);

% Store image in folder
saveas(P,[Output_path,'/Maps/Figures/Scatterer_Map/JPG/',H_Text,'.jpg'],'jpg');
saveas(P,[Output_path,'/Maps/Figures/Scatterer_Map/FIG/',H_Text,'.fig'],'fig');

%--------------------------------------------------------------------------