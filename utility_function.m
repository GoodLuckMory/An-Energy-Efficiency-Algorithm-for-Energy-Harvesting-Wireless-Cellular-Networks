function [power_sum,u_macro,u_macro_match] = utility_function(mean_power,sum_energy_per_slot,sum_energy_per_slot_updated)

time = 50;
mean_power1 = mean_power; % 高斯过程的均值
std_dev_power = sqrt(mean_power1); % 高斯过程的标准偏差

% 生成高斯过程样本
gaussian_samples = mean_power1 + std_dev_power * randn(1, time);

power_values = abs(gaussian_samples);

% 计算每个时刻的功率加上上一个时刻的功率的和
power_sum = zeros(1, time); % 初始化功率和数组
for t = 1:time
    if t == 1
        % 第一个时刻没有前一个时刻，只考虑当前时刻的功率值
        power_sum(t) = power_values(t);
    else
        % 从第二个时刻开始，功率和是当前时刻和前一个时刻的功率值之和
        power_sum(t) = power_values(t) + power_sum(t-1);
    end
end
figure;

% 目标长度
targetLength = 50;

% 计算每个分段的长度
segmentLength = length(sum_energy_per_slot) / targetLength;

% 初始化压缩后的数组
u_macro = zeros(1, targetLength);

% 对原始数组进行分段和求平均
for i = 1:targetLength
    % 计算当前段的起始和结束索引
    startIndex = round((i-1) * segmentLength) + 1;
    endIndex = round(i * segmentLength);
    
    % 计算当前段的平均值并赋值给压缩后的数组
    u_macro(i) = mean(sum_energy_per_slot(startIndex:endIndex));
end
u_macro_match = sum_energy_per_slot_updated(1:50);
end