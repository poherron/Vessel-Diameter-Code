function ZI = interp2_p(X,Y,Z,XI,YI)
%INTERP2 2-D interpolation (table lookup). CUSTOMIZED BY PYC @ MUSC.
%   ZI = INTERP2(X,Y,Z,XI,YI) interpolates to find ZI, the values of the
%   underlying 2-D function Z at the points in matrices XI and YI.
%   Matrices X and Y specify the points at which the data Z is given.
%   LINEAR INTERPOLATION, NO EXTRAPOLATION. X and Y can be VECTORS.
%   Class support for inputs X, Y, Z, XI, YI:  
%      float: double, single
%
%   See also INTERP, INTERP1, INTERP3, INTERPN, MESHGRID, GRIDDATA.
% 
% modified by Pratik Chhatbar at Kara Lab on 6/6/2012
% MUSC Neurosciences, Charleston, SC.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 5.33.4.14 $

[nrows,ncols] = size(Z);

if isempty(X), X = 1:ncols; end
if isempty(Y), Y = 1:nrows; end

s = 1 + (XI-X(1))/(X(end)-X(1))*(ncols-1);
t = 1 + (YI-Y(1))/(Y(end)-Y(1))*(nrows-1);

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if ~isempty(sout), s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if ~isempty(tout), t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
onemt = 1-t;
ZI =  ( Z(ndx).*(onemt) + Z(ndx+1).*t ).*(1-s) + ...
    ( Z(ndx+nrows).*(onemt) + Z(ndx+(nrows+1)).*t ).*s;

% Now set out of range values to NaN.
if ~isempty(sout), ZI(sout) = nan; end
if ~isempty(tout), ZI(tout) = nan; end

