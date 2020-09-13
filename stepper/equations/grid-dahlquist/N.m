function ns = N(~,yh,params)
%N Returns ODE nonlinear function
l2 = repmat(1i * params.z2(:), params.n, 1);
if(params.epsilon ~= 0)
    [l1, ~] = meshgrid(1i * params.z1(:), 1i * params.z1(:));
    l2 = l2 + params.epsilon * abs(l1(:));
end
ns = l2 .* yh;
end