function R_Map(Im_Ma,im_min,im_max,Bounce,a_p,r_p,flag_dB)

% function for displaying reflectivity maps

%--------------------------------------------------------------------------

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% input: 
% - Im_Ma: image matrix containing pixel amplitudes
% - im_min, im_max: power range in amplitude to be displayed (dB values or clipping threshold)
% - Bounce: string containing bounce level information 
% - a_p: number of pixels in azimuth
% - r_p: number of pixels in range
% - flag_dB: flag for using dB scale or clipping; 1: dB; 0: clipping

% output: none

%--------------------------------------------------------------------------

% global parameters

global r_geom range_dir ang Output_path;

% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% Output_path: absolute path to designated folder where simulation results are to be stored

%--------------------------------------------------------------------------

% Scale image entries to interval of gray values
num_g = 2^8-1; % maximum gray value --> 255
Im_Ma = im_trans(Im_Ma,im_min,im_max,num_g);

% Open new figure
figure;
set(gcf,'Name',Bounce,'Numbertitle','off')

if range_dir == 0
    Im_Ma = flipud(Im_Ma); % flip up down
    P = image(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
else
    Im_Ma = fliplr(Im_Ma); % flip left right
    P = image(a_p, r_p, Im_Ma);
    set(P,'Tag','0');
end

% Colormap --> gray values
Map = colormap(gray(num_g+1));

% Colorbar
colorbar;

if flag_dB == 1

    % display in ground range 
    if r_geom == 1
        xlabel('Azimuth Pixels'); ylabel('Ground Range Pixels'); title([Bounce ' (dB) - angle of incidence: ' num2str(ang*(180/pi)) ' degrees']);
    end

    % display in slant range
    if r_geom == 0
        xlabel('Azimuth Pixels'); ylabel('Slant Range Pixels'); title([Bounce ' (dB)']);
    end
    
else
    
    % display in ground range 
    if r_geom == 1
        xlabel('Azimuth Pixels'); ylabel('Ground Range Pixels'); title([Bounce ' (Clipping) - angle of incidence: ' num2str(ang*(180/pi)) ' degrees']);
    end

    % display in slant range
    if r_geom == 0
        xlabel('Azimuth Pixels'); ylabel('Slant Range Pixels'); title([Bounce ' (Clipping)']);
    end

end

axis equal; axis image;

% Store image in folder
saveas(P,[Output_path,'/Maps/Figures/Ref_Maps/JPG/',Bounce,'.jpg'],'jpg');
saveas(P,[Output_path,'/Maps/Figures/Ref_Maps/FIG/',Bounce,'.fig'],'fig');
imwrite(Im_Ma,Map,[Output_path,'/Maps/Frames/Ref_Maps/',Bounce,'_Fr.tif'],'Compression','none');

%--------------------------------------------------------------------------