function control = SolveLinearMPC(A_k, B_k, C_k, q, r, lower, upper, x0, refs, N_p)

    % 预测步长是N_p
    % 设状态量个数是Xn
    Xn = length(x0);
    % refs是N_p个参考状态组合成的参考状态合集
    % refs的维度是[N_p*Xn, 1]
    
    % 求矩阵k
    k = cell(N_p, N_p);
    
    for i = 1:N_p
        for j = 1:N_p
            if i < j
                k{i, j} = B_k*0;
            else
                k{i, j} = A_k^(i-j)*B_k;
            end
        end
    end
    THETA = cell2mat(k);  % 将元胞数组转换为基础数据类型的普通数组
    
    % 求矩阵M
    m = A_k*x0 + C_k;
    M = cell(N_p,1);
    
    for i = 1 : N_p
        M{i} = m;
        m = A_k * m +C_k;
    end
    
    M = cell2mat(M);
    
    % Q,R
    Q = [];
    R = [];
    
    for i = 1:N_p
        Q = blkdiag(Q, q);  % 从三个不同大小的矩阵创建一个分块对角矩阵
        R = blkdiag(R, r);  % 从三个不同大小的矩阵创建一个分块对角矩阵
    end
    
    ll = repmat(lower, N_p, 1);  % 创建一个所有元素的值均为 10 的 3×2 矩阵。A = repmat(10,3,2)
    uu = repmat(upper, N_p, 1);
    
    H = 2 * ((THETA.')*Q*THETA + R);
    f = (THETA.')*Q*(M-refs);
    
    [control,~,~,~,~] = quadprog(H, f, [],[],[],[],ll, uu);
end