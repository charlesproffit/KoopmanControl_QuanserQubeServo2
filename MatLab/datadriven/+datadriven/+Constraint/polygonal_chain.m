function constraints = polygonal_chain(Y)
%POLYGONAL_CHAIN Computes constraints for a polygonal chain in a SISO system.
% 
%   [constraints, slack] = polygonal_chain(Y)
% 
%   ----------------------------------------------------------------------------
%   Outputs
%     * constraints : List of LMI constraints in YALMIP framework.
% 
%   ----------------------------------------------------------------------------
%   References
%   [1] 
% 
%   ----------------------------------------------------------------------------
%   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
%

arguments (Input)
    Y
end

s = size(Y);
if (s(1) ~= 1) && (s(2) ~= 1)
    error("Polygonal Chain is only for SISO cases!");
end
Y = Y(:);

Y = [Y(1)'; Y; Y(end)'];
n = get_normal_direction(double(Y));

polygonalY1 = 2*real(conj(n).*Y(2:end));
polygonalY2 = 2*real(conj(n).*Y(1:end-1));

constraints = [
    polygonalY1 >= 1e-5;
    polygonalY2 >= 1e-5;
    ];

end

%% Internal Functions
function n = get_normal_direction(r)
% GET_NORMAL_DIRECTION computes the unit normal vectors of a 2D curve 
% represented as a complex-valued vector r.

arguments
    r (:, 1) double
end

% Compute initial normal directions (rotate tangent by 90° counterclockwise)
n = 1j * diff(r);

% Handle cases where n is zero
zeroIdx = (n == 0);
n(zeroIdx) = r(zeroIdx);

% Fix inconsistent normal directions
fixIdx = find(imag(conj(n).*r(1:end-1)) .* imag(conj(n).*r(2:end)) > 0);

[~, minIdx] = min(abs([ r(fixIdx), r(fixIdx+1) ]), [], 2);
n(fixIdx) = r(fixIdx - 1 + minIdx);

% Normalize the normal vectors
n = n./abs(n);

% Ensure consistency with tangent direction
n = n .* sign(real(conj(n) .* r(1:end-1)));

end