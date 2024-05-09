
%[energy_demanded_ma,energy_demanded_mi,energy_demanded_pi] = energydemand();
[power_sum,u_macro,u_macro_match,nega_energy] = algorithm(20.1,4.5,1,energy_demanded_ma,energy_demanded_mi,energy_demanded_pi);

zeta = 100; % 服务用户每耗费一焦耳能量获得的效用（单位收益）
psi = 300; % 购买一焦耳非可再生能源的价格
time = 50;
T = 1;
%figure 2
h = figure;
plot(1:time, power_sum*(zeta)*T - u_macro - nega_energy, 'r--','LineWidth', 2);
hold on
plot(1:time, power_sum*(zeta)*T - u_macro_match - nega_energy, 'r-','LineWidth', 2);
hold on
xlabel('时间');
grid on;

[power_sum,u_macro,u_macro_match,nega_energy] = algorithm(23.3,5.8,1.1,energy_demanded_ma,energy_demanded_mi,energy_demanded_pi);

figure(h)
plot(1:time, power_sum*(zeta)*T - u_macro - nega_energy , 'b--','LineWidth', 2);
hold on
plot(1:time, power_sum*(zeta)*T - u_macro_match - nega_energy , 'b-','LineWidth', 2);
hold on
legend('未利用匹配算法 第一种获取曲线','利用匹配算法 第一种获取曲线','未利用匹配算法 第二种获取曲线','利用匹配算法 第二种获取曲线','location','northwest')
xlabel('时间');
ylabel('皮蜂窝的累计效用');
grid on;