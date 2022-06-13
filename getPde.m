function aPde = getPde(p, e, t, b, f, u, time)
%getPde utility function to return or construct a pde object
%
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014 The MathWorks, Inc.

if(isa(b, 'pde.PDEModel'))
  aPde = b;
else
  N = pdeGetN( p, e, t, b, f, u, time );
  aPde = pde.PDEModel(N);
end
end