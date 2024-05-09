function [energy_demanded_ma,energy_demanded_mi,energy_demanded_pi] = energydemand()

lambda1 = 1/250000;  % 点密度1
lambda2 = 3/250000;   % 点密度2
lambda3 = 5/250000;  % 点密度3
area_size = [1350, 1350];  % 区域大小 [宽度, 高度]

[points1, num_points1] = generate_PPP(lambda1, area_size);%macro
[points2, num_points2] = generate_PPP(lambda2, area_size);%micro
[points3, num_points3] = generate_PPP(lambda3, area_size);%pico

%% macrocells
% 设置参数
num_time_slots = 600; % 时间槽的数量
mean_energy_demanded_ma = 22.5; % 能源demand的均值
std_dev_energy_demanded_ma = mean_energy_demanded_ma; % 能源收获的标准差
correlation_coefficient = 0.5; % 相关系数

% 计算基站之间的距离矩阵
distances = pdist2(points1, points1);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points1, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_demanded_ma = mean_energy_demanded_ma + std_dev_energy_demanded_ma * (L * white_noise);

%{
% 绘制能源收获数据
time = 1:num_time_slots;% 时间槽序号
figure;
for i = 1:num_points1
    plot(time, energy_demanded_ma(i,:), 'Color', rand(1,3));
    hold on;
end
xlabel('时间槽');
ylabel('能源demand量');
title('macro基站能源demand情况');
%}

%% microcells
% 设置参数
mean_energy_demanded_mi = 5; % 能源demand的均值
std_dev_energy_demanded_mi = sqrt(mean_energy_demanded_mi); % 能源收获的标准差

% 计算基站之间的距离矩阵
distances = pdist2(points2, points2);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points2, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_demanded_mi = mean_energy_demanded_mi + std_dev_energy_demanded_mi * (L * white_noise);


%% picrocells demanded
% 设置参数

mean_energy_demanded_pi = 0.85; % 能源收获的均值
std_dev_energy_demanded_pi = sqrt(mean_energy_demanded_pi); % 能源收获的标准差

%{
 生成PPP表示的基站位置
num_bs = poissrnd(lambda * num_time_slots); % 基站数量，根据PPP生成
bs_positions = rand(num_bs, 2) * 100; % 在区域[0,100]x[0,100]内生成基站的位置

%}

% 计算基站之间的距离矩阵
distances = pdist2(points3, points3);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points3, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_demanded_pi = mean_energy_demanded_pi + std_dev_energy_demanded_pi * (L * white_noise);

% 绘制能源收获数据
%{
time = 1:num_time_slots; % 时间槽序号
figure;
for i = 1:num_points3
    plot(time, energy_demanded_pi(i,:), 'Color', rand(1,3));
    hold on;
end
xlabel('时间槽');
ylabel('能源demand量');
title('picro基站能源demand情况');
%}
%end

%function harvested_energy_model(points1,points2.points3)
 
end