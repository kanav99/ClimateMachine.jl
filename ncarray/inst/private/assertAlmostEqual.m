% this function is necessary because of the limitation of matlab

function assertAlmostEqual(observed,expected)

% tolerance for testing
tol = 1e-10;

% for compatibility with matlab which does not have
% assert (OBSERVED, EXPECTED, TOL)

assert(max(abs(observed(:) - expected(:))) < tol)

% for octave prior to 3.8.0
if isempty(which('isequaln'))
  isequaln = @(x,y) isequalwithequalnans(x,y);
end

% check also NaNS
assert(isequal(isnan(observed),isnan(expected)))

end
