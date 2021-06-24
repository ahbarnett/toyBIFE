function a = logweisumexp(x,w)
% LOGWEISUMEXP   Stably compute log of a weighted sum of exponentials of a list.
%
% a = logweisumexp(x,w) returns log(sum(w.*exp(x))) computed in a stable way
%  that avoids over/underflow.
% Inputs:
%  x - row or column vector of values
%  w - vector of weights (must be same length as x)
% Output:
%  a - the log of the w-weighted sum over the elementwise exp of x
%
% When called without arguments, runs a self test.
%
% Example: see self-test code below the function body.

% Barnett 6/1/21
if nargin==0, test_logweisumexp; return; end
x = x(:); w = w(:);              % make col vecs
if numel(x)~=numel(w), error('x and w do not have the same lengths!'); end
xmax = max(x);
a = log(sum(w.*exp(x-xmax)));
a = a + xmax;

%%%%%%%%
function test_logweisumexp
x = linspace(0,10,1e3);
w = x(end:-1:1);
aex = log(sum(w.*exp(x)));   % true
xoff = -1000;
a = logweisumexp(x+xoff,w)-xoff;
relerr = abs(1-a/aex);
if relerr<1e-16*abs(xoff)      % tol needs to be more lenient ~ xoff.eps_mach
  fprintf('passed: rel err %.3g\n',relerr)
else
  fprintf('failed! rel err %.3g\n',relerr)
end
