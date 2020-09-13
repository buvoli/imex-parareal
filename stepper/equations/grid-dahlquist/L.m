function l1 = L(params)
%cs returns linear operator
[l1, ~] = meshgrid(1i * params.z1(:), 1i * params.z1(:));
l1 = l1 - params.epsilon * abs(l1); % optional repartitioning
l1 = l1(:).';
end