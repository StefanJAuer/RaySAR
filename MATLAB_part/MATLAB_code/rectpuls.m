function y = rectpuls(t,Tw)
%RECTPULS Sampled aperiodic rectangle generator.
%   RECTPULS(T) generates samples of a continuous, aperiodic,
%   unity-height rectangle at the points specified in array T, centered
%   about T=0.  By default, the rectangle has width 1.  Note that the
%   interval of non-zero amplitude is defined to be open on the right,
%   i.e., RECTPULS(-0.5)=1 while RECTPULS(0.5)=0.
%
%   RECTPULS(T,W) generates a rectangle of width W.
%
%   See also GAUSPULS, TRIPULS, PULSTRAN.

%   Author(s): D. Orofino, 4/96
%   Copyright 1988-2004 The MathWorks, Inc.
%       $Revision: 1.6.4.4 $

%error(nargchk(1,2,nargin,'struct'));
%error(narginchk(1,2));

if nargin<2, Tw=1;   end

% Returns unity in interval [-Tw/2,+Tw/2) (right side of interval is open)
% Because of numerical concerns, we use eps as the tolerance
y = abs(t)<Tw/2-eps;
y(abs(t-(-Tw/2))<eps) = 1.0;

% end of rectpuls.m
