function N = pdeGetN( p, e, t, b, f, u0, time )
%PDEGETN Examine different input pieces to find N

%       Copyright 2014 The MathWorks, Inc.

nu = length(u0);
if(isa(b, 'pde.PDEModel'))
  N = b.PDESystemSize;
elseif(nu > 1 && isnumeric(u0))
  % we must have a full u0 vector
  np = size(p, 2);
  if(rem(nu,np))
    error(message('pde:pdeGetN:invalidLengthU0'));
  end
  N = nu/np;
elseif(size(b, 1) > 1)
  % we can get N from the boundary matrix
  N = max(b(1,:));
elseif ~isempty(f)
  % number of rows of f
  f=pdetfxpd(p,t,u0,time,f);
  N = size(f,1);
elseif ~isempty(b)
  % b must be a boundary file since it isn't a boundary matrix
  % It might be a function of u but since we don't have a full u0 to pass
  % in we'll assume it is not and warn the user of this fact.
  try
    % Look at the first edge
    [q,g,h,r] = pdeexpd(p,e(:,1), u0, b);
  catch ex
    if(ischar(b))
      bcFuncName = b;
    else
      bcFuncName = func2str(b);
    end
	% We want to show the user a warning but don't want the stack trace
    warning('off', 'backtrace');
    warning(message('pde:pdeGetN:bcFuncEvalWarning', bcFuncName));
    warning('on', 'backtrace');
    rethrow(ex);
  end
  N=size(g,1);
else
  N = 1; % assume scalar case
end

end

