function ginput_Intersect(arg1,az_min,ra_min,azimuth_pixel,range_pixel,num_pix_a,num_pix_r,cube_size,trans_POV)

%GINPUT_POV Graphical input from mouse.
%   GINPUT_POV(N) gets N points from the current axes and and initiates up
%   to three plots per point by starting function Tomo.m

%   Function basec on MATLAB function GINPUT
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 5.32.4.10 $  $Date: 2007/03/03 04:45:45 $

%--------------------------------------------------------------------------

% Extended by:

% Stefan Auer 
% Remote Sensing Technology Institute
% DLR

%--------------------------------------------------------------------------

% local parameters

% Input:
% - arg1: number of samples to be taken
% - az_min: minimum in azimuth [m]
% - ra_min: minimum in slant or ground range [m]
% - azimuth_pixel: sampling stepwidth in azimuth [m]
% - range_pixel: sampling stepwidth in range [m]
% - num_pix_a: number of pixels in azimuth direction
% - num_pix_r: number of pixels in range direction
% - cube_size: size of cubes to be modeled [m]
% - trans_POV: string containing POV-Ray transformations to be compensated

% Output: None

%--------------------------------------------------------------------------

% global parameters

global range_dir r_geom ang;

% range_dir: range direction within simulated images; value 0: bottom up; value 1: top down
% r_geom: image geometry in range direction; 0: slant range, 1: ground range
% ang: incidence angle of simulated signal [rad]

%--------------------------------------------------------------------------

out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC') 
   tp = get(0,'TerminalProtocol');
else
   tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
   if nargout == 1,
      if nargin == 1,
         out1 = trmginput(arg1);
      else
         out1 = trmginput;
      end
   elseif nargout == 2 || nargout == 0,
      if nargin == 1,
         [out1,out2] = trmginput(arg1);
      else
         [out1,out2] = trmginput;
      end
      if  nargout == 0
         out1 = [ out1 out2 ];
      end
   elseif nargout == 3,
      if nargin == 1,
         [out1,out2,out3] = trmginput(arg1);
      else
         [out1,out2,out3] = trmginput;
      end
   end
else
   
   fig = gcf;
   figure(gcf);
   
   if nargin == 0
      how_many = -1;
      b = [];
   else
      how_many = arg1;
      b = [];
      if  ischar(how_many) ...
            || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
            || ~(fix(how_many) == how_many) ...
            || how_many < 0
         error('MATLAB:ginput:NeedPositiveInt', 'Requires a positive integer.')
      end
      if how_many == 0
         ptr_fig = 0;
         while(ptr_fig ~= fig)
            ptr_fig = get(0,'PointerWindow');
         end
         scrn_pt = get(0,'PointerLocation');
         loc = get(fig,'Position');
         pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
         out1 = pt(1); y = pt(2);
      elseif how_many < 0
         error('MATLAB:ginput:InvalidArgument', 'Argument must be a positive integer.')
      end
   end
   
   % Suspend figure functions
   state = uisuspend(fig);
   
   toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
   if ~isempty(toolbar)
        ptButtons = [uigettool(toolbar,'Plottools.PlottoolsOff'), ...
                     uigettool(toolbar,'Plottools.PlottoolsOn')];
        ptState = get (ptButtons,'Enable');
        set (ptButtons,'Enable','off');
   end

   set(fig,'pointer','crosshair');
   fig_units = get(fig,'units');
   char = 0;

   % We need to pump the event queue on unix
   % before calling WAITFORBUTTONPRESS 
   drawnow
   
   while how_many ~= 0
      % Use no-side effect WAITFORBUTTONPRESS
      waserr = 0;
      try
	keydown = wfbp; 
      catch
	waserr = 1;
      end
      if(waserr == 1)
         if(ishandle(fig))
            set(fig,'units',fig_units);
	    uirestore(state);
            error('MATLAB:ginput:Interrupted', 'Interrupted');
         else
            error('MATLAB:ginput:FigureDeletionPause', 'Interrupted by figure deletion');
         end
      end
      
      ptr_fig = get(0,'CurrentFigure');
      if(ptr_fig == fig)
         if keydown
            char = get(fig, 'CurrentCharacter');
            button = abs(get(fig, 'CurrentCharacter'));
            scrn_pt = get(0, 'PointerLocation');
            set(fig,'units','pixels')
            loc = get(fig, 'Position');
            % We need to compensate for an off-by-one error:
            pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
            set(fig,'CurrentPoint',pt);
         else
            button = get(fig, 'SelectionType');
            if strcmp(button,'open') 
               button = 1;
            elseif strcmp(button,'normal') 
               button = 1;
            elseif strcmp(button,'extend')
               button = 2;
            elseif strcmp(button,'alt') 
               button = 3;
            else
               error('MATLAB:ginput:InvalidSelection', 'Invalid mouse selection.')
            end
         end
         pt = get(gca, 'CurrentPoint');
         
         how_many = how_many - 1;
         
         if(char == 13) % & how_many ~= 0)
            % if the return key was pressed, char will == 13,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.  
            % If this was an early breakout, don't include
            % the <Return> key info in the return arrays.
            % We will no longer count it if it's the last input.
            break;
         end
         
         out1 = [out1;pt(1,1)];
         y = [y;pt(1,2)];
         b = [b;button];
      end
   end
   
   uirestore(state);
   if ~isempty(toolbar) && ~isempty(ptButtons)
        set (ptButtons(1),'Enable',ptState{1});
        set (ptButtons(2),'Enable',ptState{2});
   end
   set(fig,'units',fig_units);
   
   if nargout > 1
      out2 = y;
      if nargout > 2
         out3 = b;
      end
   else
      out1 = [out1 y];
   end
   
   % ----------------------------------------------------------------------
   
   % Additional part starts here!!!
   
   % Programmed by Stefan Auer, Jan. 2009
   
   % ----------------------------------------
   
   % 1.) Calculate azimuth and range interval
   
   % ----------------------------------------
            
   % A.) find center of selected pixel
   pix_az = round(out1(1));
   pix_ra = round(out1(2));
       
   % B.) Calculate corresponding coordinates in azimuth-range plane
   
   % range direction bottom up
   if range_dir == 0
       
       % Intervals
       
       % azimuth   
       az_int_min = az_min+(pix_az-1)*azimuth_pixel;
       az_int_max = az_min+(pix_az)*azimuth_pixel;
        
       % range       
       ra_int_min = ra_min + (num_pix_r-pix_ra)*range_pixel;
       ra_int_max = ra_min + (num_pix_r-pix_ra+1)*range_pixel;
       
       % change to slant range in the case of ground range
       if r_geom == 1
          ra_int_min = ra_int_min*sin(ang);
          ra_int_max = ra_int_max*sin(ang);
       end
       
   end
   
   % range direction top down
   if range_dir == 1
       
       % Intervals
       
       % azimuth        
       az_int_min = az_min+(num_pix_a-pix_az)*azimuth_pixel;
       az_int_max = az_min+(num_pix_a-pix_az+1)*azimuth_pixel;
       
       % range       
       ra_int_min = ra_min + (pix_ra-1)*range_pixel;      
       ra_int_max = ra_min + (pix_ra)*range_pixel;
       
       % change to slant range in the case of ground range
       if r_geom == 1
          ra_int_min = ra_int_min*sin(ang);
          ra_int_max = ra_int_max*sin(ang);
       end
   end
           
   % ----------------------------------------
   
   % 2.) Create Geometry File
   
   % ----------------------------------------
   
   Write_Intersect_Pts(az_int_min,az_int_max,ra_int_min,ra_int_max,cube_size,trans_POV);
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = [];

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
  h=findall(fig,'type','uimenu','accel','C');   % Disabling ^C for edit menu so the only ^C is for
  set(h,'accel','');                            % interrupting the function.
  keydown = waitforbuttonpress;
  current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
  if~isempty(current_char) && (keydown == 1)           % If the character was generated by the 
	  if(current_char == 3)                       % current keypress AND is ^C, set 'waserr'to 1
		  waserr = 1;                             % so that it errors out. 
	  end
  end
  
  set(h,'accel','C');                                 % Set back the accelerator for edit menu.
catch
  waserr = 1;
end
drawnow;
if(waserr == 1)
   set(h,'accel','C');                                % Set back the accelerator if it errored out.
   error('MATLAB:ginput:Interrupted', 'Interrupted');
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
