% u = pl(lambda)
%
% for testing, rate region is a hypersphere 
%
% lambda: dual variables
%
% Output
%  u: arc capacities
function u = pl(lambda)
if sum(lambda) == 0
  lambda = rand(size(lambda));
end
u = lambda/norm(lambda);

