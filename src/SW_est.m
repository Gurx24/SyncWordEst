clc; clear; close all;

% 同步字盲识别算法

d = 1;

%***************************** 生成数据 *****************************%

data = [];      % 帧长524bit
syn_word = double(hexToBinaryVector('FAF320', 24).');
for i = 1:800
    random_data = randi([0, 1], 500, 1);      % 生成长度为 (524-24) 的随机比特序列
    frame_data  = [syn_word; random_data];    % 将同步字和随机数据组合成完整帧
    data        = [data; frame_data];         % 组合帧生成数据
end

%***************************** 增加误码 *****************************%

ber = 0.01;                          % 设置误码率1%
em  = rand(length(data), 1) < ber;   % 根据误码率生成掩码
tb  = xor(data, em);                 % 使用 XOR 操作引入误码

%***************************** 参数设置 ****************************%

bits     = tb(20:end);                  % 截获的比特流数据
bits_len = length(bits);                % 截获比特流长度
win_len  = 2000;                        % 设置滑动窗口长度2000
win_num  = floor(bits_len/win_len);     % 向下取整得到的窗口数量
W        = zeros(win_num, win_len);     % 用来存储窗口

% 截获数据划分为窗口
for i = 0:win_num-1
    W(i+1,:) = bits(i*win_len+1:(i+1)*win_len);
end

%**************************** 识别E-SW1 *****************************%

D = cell(win_len, 1); % 用来存取窗口滑动m时同或的结果

% 将两个窗口做滑动同或运算
for m = 0:win_len-1
    W_1    = W(1, 1:win_len-m);    % 截断窗口
    W_2    = W(2, 1+m:win_len);    % 截断窗口
    D{m+1} = ~bitxor(W_1, W_2);
end

% 获取最长连1字串(无误码)
% maxCountArray = zeros(win_len, 1);      % 存储每个行向量中的最大连续1的数量
% maxPositionArray = zeros(win_len, 1);   % 存储每个行向量中最大连续1的起始位置
% for i = 1:length(D)
%     temp = D{i};        % 提取当前行向量
%     maxCount = 0;       % 初始化当前最大连续1的计数
%     currentCount = 0;   % 当前连续1的计数
%     maxPosition = 0;    % 当前最大连续1的起始位置
%     currentPosition = 0;% 当前连续1的起始位置
% 
%     % 遍历当前行向量
%     for j = 1:length(temp)
%         if temp(j) == 1
%             if currentCount == 0
%                 currentPosition = j;  % 记录当前连续1的起始位置
%             end
%             currentCount = currentCount + 1;  % 增加当前连续1的计数
%         else
%             % 检查是否需要更新最大值
%             if currentCount > maxCount
%                 maxCount = currentCount;        % 更新最大连续1的计数
%                 maxPosition = currentPosition;  % 更新最大连续1的起始位置
%             end
%             currentCount = 0;                   % 重置当前计数
%         end
%     end
% 
%     % 最后一次检查，以防数据以1结尾
%     if currentCount > maxCount
%         maxCount = currentCount;
%         maxPosition = currentPosition;
%     end
% 
%     % 存储结果
%     maxCountArray(i) = maxCount;        % 记录最大连续1的数量
%     maxPositionArray(i) = maxPosition;  % 记录最大连续1的起始位置
% end

% 获取最长公共字串
J   = 5;
th1 = ceil(J/2)+2;  % 上升沿阈值
th2 = J-1;          % 下降沿阈值
maxCountArray = zeros(win_len, 1);
maxPositionArray = zeros(win_len, 1);

for i = 1:length(D)
    temp = D{i};        % 提取当前行向量
    maxCount = 0;       % 初始化当前最大连续1的计数
    maxPosition = 0;    % 当前最大连续1的起始位置
    currentCount = 0;   % 当前连续1的计数
    currentPosition = 0;% 当前连续1的起始位置

    risePosition = length(temp)-J;% 防止未找到上升沿就找到下降沿的情况
    fallPosition = 0;

    for j=J+1 : length(temp)-J
        if sum(temp(j-J:j-1))<th1 && sum(temp(j:j+J-1))>th2
            risePosition = j;
        end
        if sum(temp(j-J+1:j))>th2 && sum(temp(j+1:j+J))<th1
            fallPosition = j;
        end
        
        if fallPosition > risePosition
            currentCount = fallPosition - risePosition + 1;

            % 检查是否需要更新最大值
            if currentCount > maxCount
                maxCount    = currentCount;  
                maxPosition = risePosition;
            end
            risePosition = length(temp)-J;
            fallPosition = 0;
            currentCount = 0;  % 重置当前计数
        end

    end

    % 最后一次检查，以防数据以1结尾
    if currentCount > maxCount
        maxCount    = currentCount;
        maxPosition = currentPosition;
    end

    % 存储结果
    maxCountArray(i)    = maxCount;        % 记录最大连续1的数量
    maxPositionArray(i) = maxPosition;     % 记录最大连续1的起始位置
end

% 获取最长公共字串的位置
[canditate_SW_len, canditate_SW_D] = max(maxCountArray);
canditate_SW_pos = maxPositionArray(canditate_SW_D);
canditate_SW = W(1, canditate_SW_pos:canditate_SW_pos+canditate_SW_len-1); % 获得候选同步字

% 绘制带有上升沿和下降沿的图像
if d>10
    data = D{97};
    time = 1:length(data);

    figure;
    stairs(time, data);
    ylim([-0.5, 1.5]);
end

% 定义扩展同步字集合
Lea        = 8;                         % 扩展长度为8bit
E_SW       = zeros(win_num, canditate_SW_len+2*Lea);
E_SW(1, :) = W(2, canditate_SW_pos+canditate_SW_D-1-Lea:canditate_SW_pos+canditate_SW_D-1+canditate_SW_len+Lea-1);
Le         = canditate_SW_len+2*Lea;    % 扩展同步字的总长度

%************************** 获取E-SW集合 ****************************%

C     = zeros(win_num, Le);     % 定义相关性矩阵
C_max = zeros(1, win_num-1);    % 各窗口对应相关系数最大值
n_pos = zeros(1, win_num-1);    % 各窗口对应相关系数最大值的位置

for i = 2:win_num
    for n = 0:win_len-Le
        C(i-1, n+1) = (1/Le)*sum(~bitxor(E_SW(1, :), W(i, n+1:n+Le)));
    end
    [C_max(i-1), n_pos(i-1)] = max(C(i-1, :));
    E_SW(i, :) = W(i, n_pos(i-1):n_pos(i-1)+Le-1);
end

% 绘制相关系数图像
% figure;
% plot(C(2, :));

%**************************** 确定SW码字内容 ************************%

garma_sum = (1/win_num)*sum(E_SW);
index     = 1;
for j = 1:Le
    if garma_sum(j) > 0.85        % 设定阈值
        est_sw(index) = 1;
        index = index + 1;
    else
        if garma_sum(j) < 0.15    % 设定阈值
            est_sw(index) = 0;
            index = index + 1;
        end
    end
end

disp(['估计出的同步字为:', binaryVectorToHex(est_sw(1:end), 'MSBFirst')]);
