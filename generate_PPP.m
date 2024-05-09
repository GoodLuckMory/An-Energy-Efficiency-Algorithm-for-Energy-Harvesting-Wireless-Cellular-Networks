function [points, num_points] = generate_PPP(lambda, area_size)
% lambda: 点密度
% area_size: 区域大小 [width, height]

% 计算区域面积
area = area_size(1) * area_size(2);

% 生成点的数量
num_points = poissrnd(lambda * area);

% 在二维空间中随机生成点的位置
points = rand(num_points, 2) .* repmat(area_size, num_points, 1);

end
