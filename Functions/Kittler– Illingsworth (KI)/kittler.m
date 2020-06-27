function [T, V] = kittler(I, h, g)
% T = KITTLER(I) computes a global threshold (T) of the image I
% according to the method described by J. KITTLER and J. ILLINGWORTH.
%
% T = KITTLER([], H, G) computes a threshold for the image with
% histogram H. The second argument (optional) G can contain graylevel 
% values for the histogram bins (defaults to 0:(length(H)-1))
%
% Author: Bob Pepin (bobpepin@gmail.com)
%
if ~isempty(I)
   h = histc(I(:), 0:0.001:ceil(max(I(:))));
end
if nargin < 2
    g = (0:(length(h)-1))';
else
    g = reshape(g, [], 1);
end
C = cumsum(h);
M = cumsum(h .* g);
S = cumsum(h .* g.^2);
sigma_f = sqrt(S ./ C - (M./C).^2);
Cb = C(end) - C;
Mb = M(end) - M;
Sb = S(end) - S;
sigma_b = sqrt(Sb ./ Cb - (Mb./Cb).^2);
P = C ./ C(end);
V = P .* log(sigma_f) + (1-P).*log(sigma_b) - P.*log(P) - (1-P).*log(1-P);
V(~isfinite(V)) = Inf;
[~, idx] = min(V);
T = g(idx)*0.001;
end