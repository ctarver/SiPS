function undock(hFig)
%UNDOCK specified figure(s).
%   UNDOCK(H) sets the WindowsStyle property of the specified figure(s) to
%       Normal. This syntax is simpler than set(h,'WindowStyle','Normal').
%
%   UNDOCK without any inputs simply undocks the current figure window.
%
%   See also DOCK.

% Copyright 2005-2010 The MathWorks, Inc.

%defensive programming
error(nargchk(0,1,nargin))
error(nargoutchk(0,0,nargout))

%default
if nargin<1
  hFig = gcf;
end

set(hFig,'WindowStyle','Normal')
