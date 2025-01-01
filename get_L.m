function [acc_l] = get_L(X, U)
N = size(X,1);
T = vec2lab(U');
acc_s = zeros(N,N);
e=7;
ker = struct('type','gauss','width',1);p=1;
parfor i = 1:N
    % 提取第 i 行数据
    x_vec = X(i,:);
    % 计算所有点到第 i 个点的距离
    x_dist = sqrt(sum((X - x_vec).^2, 2));
    % 排序并找到最近邻的索引
    [~, index] = sort(x_dist);
    j = index(1:e); % 提取 e 个最近邻的索引
    % 局部存储结果
%     if T(i,:) == T(j,:)
%         p = 1;
%     else
%         p = 1.05;
%     end
    temp_acc = zeros(1, N); % 局部变量
    for col = 1:e
        temp_acc(j(col)) = kernel(ker, X(i, :)', X(j(col), :)') * p;
    end
    % 返回结果
    result{i} = temp_acc; % 用 cell 数组存储每次迭代的结果
end
% 合并结果
parfor i = 1:N
    acc_s(i, :) = result{i};
end
acc_d = diag( sum( acc_s));
acc_l = zeros(size(acc_d));
acc_l = acc_d - acc_s;
end