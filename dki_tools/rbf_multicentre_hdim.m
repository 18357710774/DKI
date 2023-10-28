function y = rbf_multicentre_hdim(xyz, k, xc, a, delta, w)
% Y = rbf_multicentre(xyz,k,delta,xc,w)
% is a linear combination of RBF functions on sphere with centres xc and
% weigts w.
%
% Inputs:
% xyz -- points to evaluate, size(xyz,1)=Num. Points, size(xyz,2)=3
% k -- type of Wendland; smoothness = k + 3/2
% delta -- Scaling factor delta > 0 (Default delta = 1)
% xc -- set of centres of the Wendland functions;
% size(xc) = [Num. centres, Dim. sphere +1]
% w -- set of weights; row vector


if nargin < 5
    delta = 1;
end

if nargin < 6
    w = ones(size(xc, 1),1);
end

% r = xyz' * xc; % n.points X n.center
r = pdist2(xyz, xc, 'euclidean');  % n.points X n.center
r = r./a;
Y = Wendland_r(r, k, delta);
y = Y * w;
