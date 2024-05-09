function [power_sum,u_macro,u_macro_match,sum_energy_pi_per_slot50] = algorithm(mean_ma,mean_mi,mean_pi,energy_demanded_ma,energy_demanded_mi,energy_demanded_pi)

c_macro1 = 24;
c_micro1 = 7;
c_pico1 = 3;
c_macro2 = 30;
c_micro2 = 15;
c_pico2 = 8;

lambda1 = 1/250000;  % 点密度1
lambda2 = 3/250000;   % 点密度2
lambda3 = 5/250000;  % 点密度3
area_size = [1350, 1350];  % 区域大小 [宽度, 高度]
correlation_coefficient = 0.5; % 相关系数
num_time_slots = 600;

% 生成三个不同密度的PPP分布
[points1, num_points1] = generate_PPP(lambda1, area_size);%macro
[points2, num_points2] = generate_PPP(lambda2, area_size);%micro
[points3, num_points3] = generate_PPP(lambda3, area_size);%pico

%{
% 绘制图像
figure;
% 定义区域大小和透明度
area_x = [0, area_size(1), area_size(1), 0];
area_y = [0, 0, area_size(2), area_size(2)];
alpha_level = 0.2;  % 平面透明度

% 在 z = 3 的高度添加平面
fill3(area_x, area_y, [3 3 3 3], 'b', 'FaceAlpha', alpha_level, 'DisplayName', 'Plane at macrocells');
hold on;
% macrocells 在 z = 3 的高度
scatter3(points1(:, 1), points1(:, 2), 3 * ones(size(points1, 1), 1), 'b.', 'DisplayName', 'macrocells');

% 在 z = 2 的高度添加平面
fill3(area_x, area_y, [2 2 2 2], 'r', 'FaceAlpha', alpha_level, 'DisplayName', 'Plane at microcells');
% microcells 在 z = 2 的高度
scatter3(points2(:, 1), points2(:, 2), 2 * ones(size(points2, 1), 1), 'r.', 'DisplayName', 'microcells');

% 在 z = 1 的高度添加平面
fill3(area_x, area_y, [1 1 1 1], 'g', 'FaceAlpha', alpha_level, 'DisplayName', 'Plane at picocells');
% picocells 在 z = 1 的高度
scatter3(points3(:, 1), points3(:, 2), 1 * ones(size(points3, 1), 1), 'g.', 'DisplayName', 'picocells');

xlabel('X轴');
ylabel('Y轴');
zlabel('Z轴');  % 添加Z轴标签
title('不同密度的三维PPP分布');
axis([0, area_size(1), 0, area_size(2), 0, 4]);  % 添加Z轴范围
legend('Location', 'northeast');
grid on;
%}

%average users connect to bs
%avg_users_macro=floor(average_users(0.65,80/250000,40,40,6.3,1,lambda1,lambda2,lambda3,4));
%avg_users_micro=floor(average_users(0.65,80/250000,6.3,40,6.3,1,lambda1,lambda2,lambda3,4));
%avg_users_picro=floor(average_users(0.65,80/250000,1,40,6.3,1,lambda1,lambda2,lambda3,4));



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
%}
%% macrocells harvested

% 设置参数
mean_energy_harvested_ma = mean_ma; % 能源收获的均值
std_dev_energy_harvested_ma = sqrt(mean_energy_harvested_ma); % 能源收获的标准差


%{
 生成PPP表示的基站位置
num_bs = poissrnd(lambda * num_time_slots); % 基站数量，根据PPP生成
bs_positions = rand(num_bs, 2) * 100; % 在区域[0,100]x[0,100]内生成基站的位置

%}

% 计算基站之间的距离矩阵
distances = pdist2(points1, points1);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points1, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_harvested_ma = mean_energy_harvested_ma + std_dev_energy_harvested_ma * (L * white_noise);

temp = energy_harvested_ma;
%energy_harvested_ma = min(temp, c_macro1);
energy_harvested_ma2 = min(temp, c_macro2);

%% microcells
% 设置参数

mean_energy_harvested_mi = mean_mi; % 能源收获的均值
std_dev_energy_harvested_mi = sqrt(mean_energy_harvested_mi); % 能源收获的标准差

%{
 生成PPP表示的基站位置
num_bs = poissrnd(lambda * num_time_slots); % 基站数量，根据PPP生成
bs_positions = rand(num_bs, 2) * 100; % 在区域[0,100]x[0,100]内生成基站的位置

%}

% 计算基站之间的距离矩阵
distances = pdist2(points2, points2);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points2, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_harvested_mi = mean_energy_harvested_mi+ std_dev_energy_harvested_mi * (L * white_noise);
temp = energy_harvested_mi;
%energy_harvested_mi = min(temp, c_micro1);
energy_harvested_mi2 = min(temp, c_micro2);

% 绘制能源收获数据
%{
time = 1:num_time_slots; % 时间槽序号
figure;
for i = 1:num_points2
    plot(time, energy_harvested_mi(i,:), 'Color', rand(1,3));
    hold on;
end
xlabel('时间槽');
ylabel('能源收获量');
title('micro基站能源收获情况');
%}

%% picrocells
% 设置参数
mean_energy_harvested_pi = mean_pi; % 能源收获的均值
std_dev_energy_harvested_pi = sqrt(mean_energy_harvested_pi); % 能源收获的标准差

% 计算基站之间的距离矩阵
distances = pdist2(points3, points3);

% 计算基站之间的相关系数
corr_matrix = exp(-correlation_coefficient * distances);

% 进行Cholesky分解
L = chol(corr_matrix, 'lower');

% 生成随机白噪声
white_noise = randn(num_points3, num_time_slots);

% 通过Cholesky分解生成具有指定相关性的随机过程
energy_harvested_pi = mean_energy_harvested_pi + std_dev_energy_harvested_pi * (L * white_noise);
temp = energy_harvested_pi;
%energy_harvested_pi = min(temp, c_pico1);
energy_harvested_pi2 = min(temp, c_pico2);
% 绘制能源收获数据
%{
time = 1:num_time_slots; % 时间槽序号
figure;
for i = 1:num_points3
    plot(time, energy_harvested_pi(i,:), 'Color', rand(1,3));
    hold on;
end
xlabel('时间槽');
ylabel('能源收获量');
title('picro基站能源收获情况');
%}


%end

%% test

% capacity of macro
energy_surplus_ma = energy_harvested_ma - energy_demanded_ma;
temp = energy_surplus_ma;
energy_surplus_ma = min(temp,c_macro1);
energy_surplus_ma50 = abs(energy_surplus_ma(1:50));
% micro
energy_surplus_mi = energy_harvested_mi - energy_demanded_mi;
temp = energy_surplus_mi;
energy_surplus_mi = min(temp,c_micro1);
%pico
energy_surplus_pi = energy_harvested_pi - energy_demanded_pi;
temp = energy_surplus_pi;
energy_surplus_pi = min(temp,c_pico1);

%  another capacity of macro
energy_surplus_ma2 = energy_harvested_ma - energy_demanded_ma;
temp = energy_surplus_ma2;
energy_surplus_ma2 = min(temp,c_macro2);
% micro
energy_surplus_mi2 = energy_harvested_mi - energy_demanded_mi;
temp = energy_surplus_mi2;
energy_surplus_mi2 = min(temp,c_micro2);
%pico
energy_surplus_pi2 = energy_harvested_pi - energy_demanded_pi;
temp = energy_surplus_pi2;
energy_surplus_pi2 = min(temp,c_pico2);


%% macrocells buy and sell

% 初始化买方和卖方基站的数组
sellers_ma_match = {};
buyers_ma_match = {};
% 初始化买方和卖方基站编号数组，这次为每个时间槽创建一个单元数组
sellers_ma = cell(1, num_time_slots);
buyers_ma = cell(1, num_time_slots);
% 假设 num_time_slots 和 num_points1 已经定义
% 初始化累积买方能量需求总和为0
sum_energy_ma = 0; % 现在是一个标量，不再是数组

% 初始化用于存储每个时隙累积总和的数组
sum_energy_ma_per_slot = zeros(1, num_time_slots); % 每个元素代表一个时隙的累积总和

for t = 1:num_time_slots
    current_sellers = [];
    current_buyers = [];
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points1
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_ma(i, t) - energy_demanded_ma(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
            current_buyers(end+1) = i; % 添加基站编号
        elseif energy_diff > 0
            % 如果能量差为负，则为买方
            current_sellers(end+1) = i; % 添加基站编号
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_ma = sum_energy_ma + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_ma_per_slot(t) = sum_energy_ma;
    sellers_ma_match{t} = current_sellers;
    buyers_ma_match{t} = current_buyers;
end
%plot(1:num_time_slots, sum_energy_ma_per_slot, '-', 'LineWidth', 2);
%hold on 
sum_energy_ma = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points1
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_ma2(i, t) - energy_demanded_ma(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_ma = sum_energy_ma + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_ma_per_slot2(t) = sum_energy_ma;
end


%{
% 遍历每个时间槽
for t = 1:num_time_slots
    % 在当前时间槽，初始化当前买方和卖方列表
    current_sellers = [];
    current_buyers = [];
     
    % 遍历每个基站
    for i = 1:num_points1
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_ma(i, t) - energy_demanded_ma(i, t);
        
        % 根据能量差判断基站为买方还是卖方
        if energy_diff > 0
            % 如果能量差为正，则为卖方
            current_sellers(end+1) = i; % 添加基站编号
        elseif energy_diff < 0
            % 如果能量差为负，则为买方
            current_buyers(end+1) = i; % 添加基站编号
             sum_energy_ma = abs(energy_diff) + sum_energy_ma;
        end


        % 注意：如果能量差为零，这里未分类，可以根据需求调整
    end

%存储当前时间槽的买方和卖方基站编号
    sellers_ma{t} = current_sellers;
    buyers_ma{t} = current_buyers;
end 
%}
%{
% 显示买方和卖方基站的编号
% 由于现在有多个时间槽，您可以选择查看特定时间槽或遍历所有时间槽
for t = 1:num_time_slots
    fprintf('时间槽 %d 的卖方基站编号：\n', t);
    disp(sellers{t});
    
    fprintf('时间槽 %d 的买方基站编号：\n', t);
    disp(buyers{t});
end
%}

%% microcells buy and sell


% 初始化买方和卖方基站的数组
sellers_mi_match = {};
buyers_mi_match = {};
% 初始化买方和卖方基站编号数组，这次为每个时间槽创建一个单元数组
sellers_mi = cell(1, num_time_slots);
buyers_mi = cell(1, num_time_slots);
% 假设 num_time_slots 和 num_points1 已经定义
% 初始化累积买方能量需求总和为0
sum_energy_mi = 0; % 现在是一个标量，不再是数组

% 初始化用于存储每个时隙累积总和的数组
sum_energy_mi_per_slot = zeros(1, num_time_slots); % 每个元素代表一个时隙的累积总和

for t = 1:num_time_slots
    current_sellers = [];
    current_buyers = [];
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points2
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_mi(i, t) - energy_demanded_mi(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
            current_buyers(end+1) = i; % 添加基站编号
        elseif energy_diff > 0
            % 如果能量差为负，则为买方
            current_sellers(end+1) = i; % 添加基站编号
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_mi = sum_energy_mi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_mi_per_slot(t) = sum_energy_mi;
    sellers_mi_match{t} = current_sellers;
    buyers_mi_match{t} = current_buyers;
end
%plot(1:num_time_slots, sum_energy_ma_per_slot, '-', 'LineWidth', 2);
%hold on 
sum_energy_mi = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points2
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_mi2(i, t) - energy_demanded_mi(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_mi = sum_energy_mi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_mi_per_slot2(t) = sum_energy_mi;
end

%% picocells buy and sell


%{
% 初始化买方和卖方基站的数组
sellers_pi_match = {};
buyers_pi_match = {};

% 初始化买方和卖方基站编号数组，这次为每个时间槽创建一个单元数组
sellers_pi = cell(1, num_time_slots);
buyers_pi = cell(1, num_time_slots);
sum_energy_pi = 0;
% 遍历每个时间槽
for t = 1:num_time_slots
    % 在当前时间槽，初始化当前买方和卖方列表
    current_sellers = [];
    current_buyers = [];
   
    % 遍历每个基站
    for i = 1:num_points3
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_pi(i, t) - energy_demanded_pi(i, t);
        
        % 根据能量差判断基站为买方还是卖方
        if energy_diff > 0
            % 如果能量差为正，则为卖方
            current_sellers(end+1) = i; % 添加基站编号
        elseif energy_diff < 0
            % 如果能量差为负，则为买方
            current_buyers(end+1) = i; % 添加基站编号
            sum_energy_pi = abs(energy_diff) + sum_energy_pi;
        end
        % 注意：如果能量差为零，这里未分类，可以根据需求调整
    end
    
    % 存储当前时间槽的买方和卖方基站编号
    sellers_pi_match{t} = current_sellers;
    buyers_pi_match{t} = current_buyers;
end
%}

% 初始化买方和卖方基站的数组
sellers_pi_match = {};
buyers_pi_match = {};
% 初始化买方和卖方基站编号数组，这次为每个时间槽创建一个单元数组
sellers_pi = cell(1, num_time_slots);
buyers_pi = cell(1, num_time_slots);
sum_energy_pi = 0; % 现在是一个标量，不再是数组

% 初始化用于存储每个时隙累积总和的数组
sum_energy_pi_per_slot = zeros(1, num_time_slots); % 每个元素代表一个时隙的累积总和

for t = 1:num_time_slots
    current_sellers = [];
    current_buyers = [];
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points3
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_pi(i, t) - energy_demanded_pi(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
            current_buyers(end+1) = i; % 添加基站编号
        elseif energy_diff > 0
            % 如果能量差为负，则为买方
            current_sellers(end+1) = i; % 添加基站编号
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_pi
    sum_energy_pi = sum_energy_pi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_pi_per_slot(t) = sum_energy_pi;
    sellers_pi_match{t} = current_sellers;
    buyers_pi_match{t} = current_buyers;
end
sum_energy_pi = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points3
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_pi2(i, t) - energy_demanded_pi(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_pi
    sum_energy_pi = sum_energy_pi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_pi_per_slot2(t) = sum_energy_pi;
end

%% test

%{
% 匹配后
% 遍历每个时间点
for t = 1:num_time_slots
     % 从 cell 数组中提取当前时间点的买家索引
    current_buyers = cell2mat(buyers_ma_match(:, t));
    current_sellers = cell2mat(sellers_ma_match(:, t));
    for i = 1:length(current_buyers)
        b = current_buyers(i); % 现在 b 是一个数值，可以用作索引
        needed_energy = -energy_surplus_ma(b, t); % 买家所需能量
        % 遍历卖家直到满足当前买家需求
        for k = 1:length(current_sellers)
            s = current_sellers(k);
            if energy_surplus_ma(s, t) > 0  % 确保卖家有剩余能量
                available_energy = min(energy_surplus_ma(s, t), needed_energy);
                
                % 进行能量交换
                energy_surplus_ma(b, t) = energy_surplus_ma(b, t) + available_energy;
                energy_surplus_ma(s, t) = energy_surplus_ma(s, t) - available_energy;
                needed_energy = needed_energy - available_energy;  % 更新买家所需能量
                
                % 检查买家需求是否已满足
                if needed_energy <= 0
                    break;  % 跳出卖家循环
                end
            end
        end
    end
end

% 更新后的energy_surplus_ma反映交易后状态
% 遍历每个时间点
for t = 1:num_time_slots
     % 从 cell 数组中提取当前时间点的买家索引
    current_buyers = cell2mat(buyers_ma_match(:, t));
    current_sellers = cell2mat(sellers_ma_match(:, t));
    for i = 1:length(current_buyers)
        b = current_buyers(i); % 现在 b 是一个数值，可以用作索引
        needed_energy = -energy_surplus_ma2(b, t); % 买家所需能量
        % 遍历卖家直到满足当前买家需求
        for k = 1:length(current_sellers)
            s = current_sellers(k);
            if energy_surplus_ma2(s, t) > 0  % 确保卖家有剩余能量
                available_energy = min(energy_surplus_ma2(s, t), needed_energy);
                
                % 进行能量交换
                energy_surplus_ma2(b, t) = energy_surplus_ma2(b, t) + available_energy;
                energy_surplus_ma2(s, t) = energy_surplus_ma2(s, t) - available_energy;
                needed_energy = needed_energy - available_energy;  % 更新买家所需能量
                
                % 检查买家需求是否已满足
                if needed_energy <= 0
                    break;  % 跳出卖家循环
                end
            end
        end
    end
end

% 假设 energy_surplus_ma 是一个 600*n 的矩阵，其中行代表时间点
current_slot_sum = 0;
sum_energy_ma_test = 0;
% 遍历每个时间点
for t = 1:num_time_slots
    % 计算当前时间点所有负数值的和
    for i = 1:num_points1
        if energy_surplus_ma(i,t) < 0
            current_slot_sum = energy_surplus_ma(i,t) + current_slot_sum;
        % 累积和更新：当前时间点的负数和加上之前时间点的累积和
        end
    
    end
    sum_energy_ma_test = sum_energy_ma_test + current_slot_sum;
    sum_energy_ma_per_slot_updated_test(t) = -sum_energy_ma_test;
end

%micro test
% 匹配后
% 遍历每个时间点
for t = 1:num_time_slots
     % 从 cell 数组中提取当前时间点的买家索引
    current_buyers = cell2mat(buyers_mi_match(:, t));
    current_sellers = cell2mat(sellers_mi_match(:, t));
    for i = 1:length(current_buyers)
        b = current_buyers(i); % 现在 b 是一个数值，可以用作索引
        needed_energy = -energy_surplus_mi(b, t); % 买家所需能量
        % 遍历卖家直到满足当前买家需求
        for k = 1:length(current_sellers)
            s = current_sellers(k);
            if energy_surplus_mi(s, t) > 0  % 确保卖家有剩余能量
                available_energy = min(energy_surplus_mi(s, t), needed_energy);
                
                % 进行能量交换
                energy_surplus_mi(b, t) = energy_surplus_mi(b, t) + available_energy;
                energy_surplus_mi(s, t) = energy_surplus_mi(s, t) - available_energy;
                needed_energy = needed_energy - available_energy;  % 更新买家所需能量
                
                % 检查买家需求是否已满足
                if needed_energy <= 0
                    break;  % 跳出卖家循环
                end
            end
        end
    end
end

% 更新后的energy_surplus_ma反映交易后状态
% 遍历每个时间点
for t = 1:num_time_slots
     % 从 cell 数组中提取当前时间点的买家索引
    current_buyers = cell2mat(buyers_mi_match(:, t));
    current_sellers = cell2mat(sellers_mi_match(:, t));
    for i = 1:length(current_buyers)
        b = current_buyers(i); % 现在 b 是一个数值，可以用作索引
        needed_energy = -energy_surplus_mi2(b, t); % 买家所需能量
        % 遍历卖家直到满足当前买家需求
        for k = 1:length(current_sellers)
            s = current_sellers(k);
            if energy_surplus_mi2(s, t) > 0  % 确保卖家有剩余能量
                available_energy = min(energy_surplus_mi2(s, t), needed_energy);
                
                % 进行能量交换
                energy_surplus_mi2(b, t) = energy_surplus_mi2(b, t) + available_energy;
                energy_surplus_mi2(s, t) = energy_surplus_mi2(s, t) - available_energy;
                needed_energy = needed_energy - available_energy;  % 更新买家所需能量
                
                % 检查买家需求是否已满足
                if needed_energy <= 0
                    break;  % 跳出卖家循环
                end
            end
        end
    end
end

% 假设 energy_surplus_ma 是一个 600*n 的矩阵，其中行代表时间点
current_slot_sum = 0;
sum_energy_mi_test = 0;
% 遍历每个时间点
for t = 1:num_time_slots
    % 计算当前时间点所有负数值的和
    for i = 1:num_points2
        if energy_surplus_mi(i,t) < 0
            current_slot_sum = energy_surplus_mi(i,t) + current_slot_sum;
        % 累积和更新：当前时间点的负数和加上之前时间点的累积和
        end
    
    end
    sum_energy_mi_test = sum_energy_mi_test + current_slot_sum;
    sum_energy_mi_per_slot_updated_test(t) = -sum_energy_mi_test;
end

sum_energy_per_slot_updated_test = sum_energy_ma_per_slot_updated_test + sum_energy_mi_per_slot_updated_test;
%sum_energy_per_slot_updated_test2;
% 绘制累积和随时间变化的图表
plot(1:num_time_slots,sum_energy_per_slot_updated_test,'-');
hold on 
%plot(1:num_time_slots,sum_energy_per_slot_updated_test2,'-');
title('Cumulative Sum of Negative Values Over Time');
xlabel('Time');
ylabel('Cumulative Negative Sum');
grid on; % 添加网格线以便更清晰地查看图表细节
%}

%% asssume macrocells matching
% 初始化匹配结果存储，每个时间槽一个单元数组
matches_ma = cell(1, num_time_slots);

% 遍历每个时间槽进行匹配
for t = 1:num_time_slots
    % 当前时间槽的匹配结果
    current_matches = [];
    
    % 获取当前时间槽的买方和卖方编号
    current_sellers = sellers_ma_match{t};
    current_buyers = buyers_ma_match{t};
    
    % 为每个买方寻找最匹配的卖方
    for buyer_index = current_buyers
        min_diff = inf; % 初始化最小差值为无穷大
        match_seller = 0; % 初始化匹配的卖方编号
        
        % 获取买方的需求能量
        buyer_demand = -energy_harvested_pi(buyer_index, t) + energy_demanded_pi(buyer_index, t);
        
        % 遍历卖方寻找最佳匹配
        for seller_index = current_sellers
            % 获取卖方的供应能量
            seller_supply = energy_harvested_pi(seller_index, t) - energy_demanded_pi(seller_index, t);
            
            % 计算和买方的能量差
            diff = abs(seller_supply - buyer_demand);
            
            % 更新最小差值和匹配的卖方编号
            if diff < min_diff
                min_diff = diff;
                match_seller = seller_index;
            end
        end
        
        % 如果找到匹配的卖方，记录匹配结果
        if match_seller > 0
            current_matches = [current_matches; [buyer_index, match_seller]];
            % 一旦匹配，将卖方从当前卖方列表中移除，避免重复匹配
            current_sellers(current_sellers == match_seller) = [];
        end
    end
    
    % 存储当前时间槽的匹配结果
    matches_ma{t} = current_matches;
end

%{
for t = 1:num_time_slots
    fprintf('时间槽 %d 的匹配结果（买方基站，卖方基站）:\n', t);
    disp(matches_ma{t});
end
%}

%% asssume microcells matching
% 初始化匹配结果存储，每个时间槽一个单元数组
matches_mi = cell(1, num_time_slots);

% 遍历每个时间槽进行匹配
for t = 1:num_time_slots
    % 当前时间槽的匹配结果
    current_matches = [];
    
    % 获取当前时间槽的买方和卖方编号
    current_sellers = sellers_mi_match{t};
    current_buyers = buyers_mi_match{t};
    
    % 为每个买方寻找最匹配的卖方
    for buyer_index = current_buyers
        min_diff = inf; % 初始化最小差值为无穷大
        match_seller = 0; % 初始化匹配的卖方编号
        
        % 获取买方的需求能量
        buyer_demand = -energy_harvested_mi(buyer_index, t) + energy_demanded_mi(buyer_index, t);
        
        % 遍历卖方寻找最佳匹配
        for seller_index = current_sellers
            % 获取卖方的供应能量
            seller_supply = energy_harvested_mi(seller_index, t) - energy_demanded_mi(seller_index, t);
            % 计算和买方的能量差
            diff = abs(seller_supply - buyer_demand);
            
            % 更新最小差值和匹配的卖方编号
            if diff < min_diff
                min_diff = diff;
                match_seller = seller_index;
            end
        end
        
        % 如果找到匹配的卖方，记录匹配结果
        if match_seller > 0
            current_matches = [current_matches; [buyer_index, match_seller]];
            % 一旦匹配，将卖方从当前卖方列表中移除，避免重复匹配
            current_sellers(current_sellers == match_seller) = [];
        end
    end
    
    % 存储当前时间槽的匹配结果
    matches_mi{t} = current_matches;
end

%{
for t = 1:num_time_slots
    fprintf('时间槽 %d 的匹配结果（买方基站，卖方基站）:\n', t);
    disp(matches_mi{t});
end
%}


%% asssume picocells matching
% 初始化匹配结果存储，每个时间槽一个单元数组
matches_pi = cell(1, num_time_slots);

% 遍历每个时间槽进行匹配
for t = 1:num_time_slots
    % 当前时间槽的匹配结果
    current_matches = [];
    
    % 获取当前时间槽的买方和卖方编号
    current_sellers = sellers_pi_match{t};
    current_buyers = buyers_pi_match{t};
    
    % 为每个买方寻找最匹配的卖方
    for buyer_index = current_buyers
        min_diff = inf; % 初始化最小差值为无穷大
        match_seller = 0; % 初始化匹配的卖方编号
        
        % 获取买方的需求能量
        buyer_demand = -energy_harvested_pi(buyer_index, t) + energy_demanded_pi(buyer_index, t);
        
        % 遍历卖方寻找最佳匹配
        for seller_index = current_sellers
            % 获取卖方的供应能量
            seller_supply = energy_harvested_pi(seller_index, t) - energy_demanded_pi(seller_index, t);
            
            % 计算和买方的能量差
            diff = abs(seller_supply - buyer_demand);
            
            % 更新最小差值和匹配的卖方编号
            if diff < min_diff
                min_diff = diff;
                match_seller = seller_index;
            end
        end
        
        % 如果找到匹配的卖方，记录匹配结果
        if match_seller > 0
            current_matches = [current_matches; [buyer_index, match_seller]];
            % 一旦匹配，将卖方从当前卖方列表中移除，避免重复匹配
            current_sellers(current_sellers == match_seller) = [];
        end
    end
    
    % 存储当前时间槽的匹配结果
    matches_pi{t} = current_matches;
end

%{
for t = 1:num_time_slots
    fprintf('时间槽 %d 的匹配结果（买方基站，卖方基站）:\n', t);
    disp(matches_pi{t});
end
%}

%% trade macro
energy_diff_ma = energy_harvested_ma - energy_demanded_ma;

% 对每个时间槽进行处理
for t = 1:num_time_slots
    % 获取当前时间槽的所有匹配对
    current_matches = matches_ma{t};
    
    % 遍历当前时间槽的所有匹配进行交易
    for i = 1:size(current_matches, 1)
        buyer_idx = current_matches(i, 1); % 买方索引
        seller_idx = current_matches(i, 2); % 卖方索引


        % 交易量为买方需求与卖方供应的最小值
        trade_amount = min(-energy_diff_ma(buyer_idx, t), energy_diff_ma(seller_idx, t));

        % 更新买卖双方的能量差
        energy_diff_ma(buyer_idx, t) = energy_diff_ma(buyer_idx, t) + trade_amount; % 买方减少需求
        energy_diff_ma(seller_idx, t) = energy_diff_ma(seller_idx, t) - trade_amount; % 卖方减少供应
    end
end

% 更新能量状态表格（例如，对于需求能量表格）
% 如果 energy_diff_ma 为负数，则需求未完全满足，需要从不可再生资源中获取
energy_demanded_ma_updated = energy_demanded_ma + min(energy_diff_ma, 0); 

sum_energy_ma = 0;

%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points1
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_ma(i, t) - energy_demanded_ma_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_ma = sum_energy_ma + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_ma_per_slot_updated(t) = sum_energy_ma;
end
sum_energy_ma = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points1
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_ma2(i, t) - energy_demanded_ma_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_ma
    sum_energy_ma = sum_energy_ma + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_ma_per_slot_updated2(t) = sum_energy_ma;
end
%% trade micro
energy_diff_mi = energy_harvested_mi - energy_demanded_mi;

% 对每个时间槽进行处理
for t = 1:num_time_slots
    % 获取当前时间槽的所有匹配对
    current_matches = matches_ma{t};
    
    % 遍历当前时间槽的所有匹配进行交易
    for i = 1:size(current_matches, 1)
        buyer_idx = current_matches(i, 1); % 买方索引
        seller_idx = current_matches(i, 2); % 卖方索引


        % 交易量为买方需求与卖方供应的最小值
        trade_amount = min(-energy_diff_mi(buyer_idx, t), energy_diff_mi(seller_idx, t));

        % 更新买卖双方的能量差
        energy_diff_mi(buyer_idx, t) = energy_diff_mi(buyer_idx, t) + trade_amount; % 买方减少需求
        energy_diff_mi(seller_idx, t) = energy_diff_mi(seller_idx, t) - trade_amount; % 卖方减少供应
    end
end

% 更新能量状态表格（例如，对于需求能量表格）
% 如果 energy_diff_ma 为负数，则需求未完全满足，需要从不可再生资源中获取
energy_demanded_mi_updated = energy_demanded_mi + min(energy_diff_mi, 0); 

sum_energy_mi = 0;

%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points2
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_mi(i, t) - energy_demanded_mi_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_mi
    sum_energy_mi = sum_energy_mi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_mi_per_slot_updated(t) = sum_energy_mi;
end
sum_energy_mi = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points2
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_mi2(i, t) - energy_demanded_mi_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_mi
    sum_energy_mi = sum_energy_mi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_mi_per_slot_updated2(t) = sum_energy_mi;
end



%% trade picro
energy_diff_pi = energy_harvested_pi - energy_demanded_pi;

% 对每个时间槽进行处理
for t = 1:num_time_slots
    % 获取当前时间槽的所有匹配对
    current_matches = matches_pi{t};
    
    % 遍历当前时间槽的所有匹配进行交易
    for i = 1:size(current_matches, 1)
        buyer_idx = current_matches(i, 1); % 买方索引
        seller_idx = current_matches(i, 2); % 卖方索引


        % 交易量为买方需求与卖方供应的最小值
        trade_amount = min(-energy_diff_pi(buyer_idx, t), energy_diff_pi(seller_idx, t));

        % 更新买卖双方的能量差
        energy_diff_pi(buyer_idx, t) = energy_diff_pi(buyer_idx, t) + trade_amount; % 买方减少需求
        energy_diff_pi(seller_idx, t) = energy_diff_pi(seller_idx, t) - trade_amount; % 卖方减少供应
    end
end

% 更新能量状态表格（例如，对于需求能量表格）
% 如果 energy_diff_ma 为负数，则需求未完全满足，需要从不可再生资源中获取
energy_demanded_pi_updated = energy_demanded_pi + min(energy_diff_pi, 0); 

sum_energy_pi = 0;

%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points3
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_pi(i, t) - energy_demanded_pi_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_pi
    sum_energy_pi = sum_energy_pi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_pi_per_slot_updated(t) = sum_energy_pi;
end
sum_energy_pi = 0;
%不同电池容量
for t = 1:num_time_slots
    current_slot_sum = 0; % 初始化当前时隙的买方能量需求总和为0
    % 遍历每个基站
    for i = 1:num_points3
        % 计算第i个基站在第t个时间槽的能量差
        energy_diff = energy_harvested_pi2(i, t) - energy_demanded_pi_updated(i, t);
        
        % 如果能量差为负，则为买方，计算其能量需求
        if energy_diff < 0
            current_slot_sum = current_slot_sum + abs(energy_diff);
        end
    end
    % 累加当前时隙的买方能量需求总和到sum_energy_pi
    sum_energy_pi = sum_energy_pi + current_slot_sum;
    
    % 存储累积总和到对应的时隙中
    sum_energy_pi_per_slot_updated2(t) = sum_energy_pi;
end

sum_energy_per_slot = sum_energy_ma_per_slot + sum_energy_mi_per_slot + sum_energy_pi_per_slot;
sum_energy_per_slot2 = sum_energy_ma_per_slot2;
sum_energy_per_slot_updated(t) = sum_energy_ma_per_slot_updated(t);
sum_energy_per_slot_updated2(t) = sum_energy_ma_per_slot_updated2(t)...
+sum_energy_mi_per_slot_updated2(t)+sum_energy_pi_per_slot_updated2(t);



%% utility function_

[power_sum,u_macro,u_macro_match] = utility_function(mean_pi,sum_energy_pi_per_slot,sum_energy_pi_per_slot_updated);
sum_energy_ma_per_slot50 = sum_energy_ma_per_slot(1:50);
sum_energy_pi_per_slot50 = sum_energy_pi_per_slot(1:50);

%% overall
% 绘制每个时隙的累积买方能量需求总和随时间增长的图像
% figure 1
plot(1:num_time_slots, sum_energy_per_slot, 'b-', 'LineWidth', 2);
hold on
plot(1:num_time_slots, sum_energy_per_slot2, 'r--', 'LineWidth', 2);
hold on
plot(1:num_time_slots, sum_energy_ma_per_slot_updated, 'g-', 'LineWidth', 2);
hold on
plot(1:num_time_slots, sum_energy_pi_per_slot_updated2, 'black--', 'LineWidth', 2);
hold on
xlabel('时间');
ylabel('已使用的不可再生资源总量');
legend('未使用匹配算法 不同电池容量c1=24, c2=7, c3=3','未使用匹配算法 不同电池容量c1=30, c2=10, c3=8','使用匹配算法 不同电池容量c1=24,c2=7 c3=3','使用匹配算法 不同电池容量c1=30,c2=10,c3=8','location','northwest');
grid on;

figure;
plot(1:num_time_slots, -(sum_energy_per_slot_updated - sum_energy_per_slot), 'g-', 'LineWidth', 2);
hold on
plot(1:num_time_slots, -(sum_energy_per_slot_updated2 - sum_energy_per_slot2), 'yellow--', 'LineWidth', 2);
hold on
xlabel('时间');
ylabel('能源效率提升量');
legend('不同电池容量c1=24, c2=7, c3=3','不同电池容量c1=30, c2=10, c3=8','location','northwest');
grid on;




save para;

end
