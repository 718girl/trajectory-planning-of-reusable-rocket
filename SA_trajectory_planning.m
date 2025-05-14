function [state_history, control_history] = SA_rocket_trajectory()
    tic;
    %% 参数初始化
    N = 50;            % 时间节点数
    dt = 0.8;            % 时间步长(s)
    T0 = 1000;         % 初始温度
    alpha = 0.95;      % 降温系数
    max_iter = 50;     % 每温度迭代次数
    
    %% 控制变量边界 (dλ/dt, dμ1/dt, dμ2/dt, dφ/dt, dψ/dt)
    control_lim = [-0.3,  0.3;       % dλ/dt 
                  -5, 5;      % dμ1/dt [°/s]
                  -5, 5;      % dμ2/dt [°/s]
                  -5, 5;    % dφ/dt [°/s]
                  -5, 5]';  % dψ/dt [°/s]
    
    %% 模拟退火初始化
    current_controls = initialize_controls(N, control_lim);
    [current_cost, current_state] = evaluate_trajectory(current_controls, dt);
    
    %% 退火主循环
    T = T0;
    while T > 1e-3
        for iter = 1:max_iter
            % 生成邻域解
            new_controls = generate_neighbor(current_controls, control_lim, T/T0);
            
            % 评估轨迹
            [new_cost, new_state] = evaluate_trajectory(new_controls, dt);
            
            % Metropolis准则
            if accept_solution(current_cost, new_cost, T)
                current_controls = new_controls;
                current_cost = new_cost;
                current_state = new_state;
            end
        end
        T = T * alpha;
    end
    
    %% 输出结果
    total_time = toc;
    fprintf('总运行时间为%.5f\n', total_time);
    [cost, state_terminal, state_history, control_history] = evaluate_trajectory(current_controls, dt);
    fprintf('代价函数结果%.5f\n', cost);
    save_trajectory_data(state_history, control_history, dt)
    visualize_results(state_history, control_history, dt);
end

%% 关键子函数
function controls = initialize_controls(N, lim)
    dim = size(lim,2);
    controls = zeros(dim, N);
    for i = 1:dim
        controls(i,:) = lim(1,i) + (lim(2,i)-lim(1,i))*rand(1,N);
    end
end

function new_controls = generate_neighbor(old_controls, lim, temp_factor)
    perturbation = diag(lim(2,:)-lim(1,:)) * 0.1 * temp_factor;
    new_controls = old_controls + perturbation*randn(size(old_controls));
    
    % 施加边界约束
    for i = 1:size(lim,2)
        new_controls(i,:) = max(lim(1,i), min(lim(2,i), new_controls(i,:)));
    end
end

function [cost, state, state_history, control_history] = evaluate_trajectory(controls, dt)
    % 状态变量初始化
    state.rx = -600;        % X位置(m)
    state.ry = 6000;        % Y高度(m)
    state.rz = -200;        % Z位置(m)
    state.vx = -80;         % X速度(m/s)
    state.vy = -300;        % Y速度(m/s)
    state.vz = 30;          % Z速度(m/s)
    state.m = 49000;        % 质量(kg)
    state.lambda = 0.35;    % 发动机开度（初始化为下限）
    state.mu1 = 5;          % 发动机摆角(°)
    state.mu2 = 5;          % 发动机摆角(°)
    state.phi = 80;         % 俯仰角(°)
    state.psi = -3;         % 偏航角(°)
    
    n_steps = size(controls,2);
    state_history = repmat(state, [1 n_steps]);
    control_history = zeros(5, n_steps);
    penalty = 0; % 过程约束惩罚

    for k = 1:n_steps
        % ================== 过程约束处理 ==================
        % 约束1: 发动机开度lambda ∈ [0.35,1]（硬约束）
        state.lambda = max(0.35, min(1, state.lambda));
        
        % 约束2: 发动机摆角mu1,mu2 ∈ [-5,5]°（硬约束）
        state.mu1 = max(-5, min(5, state.mu1));
        state.mu2 = max(-5, min(5, state.mu2));
        
        % 约束3: phi ∈ [60,90]°（硬约束）
        state.phi = max(60, min(90, state.phi));
        
        % 约束4: psi ∈ [-5,5]°（硬约束）
        state.psi = max(-5, min(5, state.psi));
        
        % 约束5: 质量m ≥ 28000kg（硬约束）
        if state.m < 28000
            cost = inf;
            return; % 立即终止计算
        end
        
        % 约束6: 高度单调递减（软约束）
        if k > 1 && state.ry > state_history(k-1).ry
            penalty = penalty + 1e5 * (state.ry - state_history(k-1).ry);
        end
        
        % ================== 控制量更新 ==================
        raw_dlambda = controls(1,k);
        raw_dmu1 = controls(2,k);
        raw_dmu2 = controls(3,k);
        raw_dphi = controls(4,k);
        raw_dpsi = controls(5,k);
        
        % ================== 过程约束处理 ==================
        % 约束处理后的控制量
        [dlambda, dmu1, dmu2, dphi, dpsi] = apply_constraints(...
            raw_dlambda, raw_dmu1, raw_dmu2, raw_dphi, raw_dpsi, state);
        
        % 记录实际应用的控制量
        control_history(:,k) = [dlambda; dmu1; dmu2; dphi; dpsi]; % 存储约束后的控制量
        
        % 更新控制相关状态（已施加硬约束）
        state.lambda = state.lambda + dlambda*dt;
        state.mu1 = state.mu1 + dmu1*dt;
        state.mu2 = state.mu2 + dmu2*dt;
        state.phi = state.phi + dphi*dt;
        state.psi = state.psi + dpsi*dt;
        
        % ================== 动力学计算 ==================
        Tmax = 912000 - 67000 * exp(-state.ry/7110);
        T = state.lambda * Tmax;
        
        % 推力分量计算（使用角度约束后的值）
        Tx_body = T * cosd(state.mu1) * cosd(state.mu2);
        Ty_body = T * cosd(state.mu1) * sind(state.mu2);
        Tz_body = -T * sind(state.mu1);
        
%         Fx = Tx_body * cosd(state.phi)*cosd(state.psi) - Ty_body*sind(state.phi)*cosd(state.psi) + Tz_body*sind(state.psi);
%         Fy = Tx_body * sind(state.phi) + Ty_body * cosd(state.phi);
%         Fz = -Tx_body * cosd(state.phi)*sind(state.psi) + Ty_body*sind(state.phi)*sind(state.psi) - Tz_body*cosd(state.psi);
        Fx = Tx_body .* cosd(state.phi) .* cosd(state.psi) - Ty_body .* sind(state.phi) + Tz_body .*cosd(state.phi).* sind(state.psi);
        Fy = Tx_body .* sind(state.phi) .* cosd(state.psi) + Ty_body .* cosd(state.phi) + Tz_body .*sind(state.phi).* sind(state.psi);
        Fz = -Tx_body .* sind(state.psi) + Tz_body .* cosd(state.psi);
        
        rho = 1.225 * exp(-state.ry/7110);
        V = norm([state.vx, state.vy, state.vz]);
        Dx = -0.5 * rho * 10.52 * 1.5 * V * state.vx;
        Dy = -0.5 * rho * 10.52 * 1.5 * V * state.vy;
        Dz = -0.5 * rho * 10.52 * 1.5 * V * state.vz;
        
        isp = 312 - 29 * exp(-state.ry/7110);
        
        % 状态更新
        state.rx = state.rx + state.vx*dt;
        state.ry = state.ry + state.vy*dt;
        state.rz = state.rz + state.vz*dt;
        state.vx = state.vx + (Fx + Dx)/state.m * dt;
        state.vy = state.vy + (Fy + Dy)/state.m * dt - 9.807*dt;
        state.vz = state.vz + (Fz + Dz)/state.m * dt;
        state.m = state.m - T/(isp*9.807)*dt;
        
        state_history(k) = state;
    end
    
    % ================== 终端约束处理 ==================
    terminal_penalty = 0;
    
    % 高度约束 ry=18.6±0.1m
    if abs(state.ry - 18.6) > 0.1
        terminal_penalty = terminal_penalty + 1e8*(state.ry - 18.6)^2; %8 
    end
    
    % 位置约束 rx,rz ∈ [-25,25]m
    if abs(state.rx) > 25 || abs(state.rz) >25
        terminal_penalty = terminal_penalty + 1e8*(max(abs(state.rx),abs(state.rz))-25)^2; %7 保证终端位置8
    end
    
    % 姿态约束 phi=90±1°, psi=0±0.5°
    if abs(state.phi - 90) > 1
        terminal_penalty = terminal_penalty + 1e10*(state.phi - 90)^2;
    end
    if abs(state.psi) > 0.5
        terminal_penalty = terminal_penalty + 1e10*(state.psi)^2;
    end
    
    % 速度约束 vx,vz∈[-3,3], vy∈[-1,1] m/s
    if abs(state.vx) >3 || abs(state.vz)>3
        terminal_penalty = terminal_penalty + 1e9*(max(abs(state.vx),abs(state.vz))-3)^2;  %10 保证速度11
    end
    if abs(state.vy) >1
        terminal_penalty = terminal_penalty + 1e9*(abs(state.vy)-1)^2;
    end
    
    % 综合代价计算
    cost = -state.m + penalty + terminal_penalty;
    
    % 可视化检查（调试时启用）
    if rand < 0.01 % 1%概率抽样检查
        fprintf('终端状态检查: rx=%.2f, rz=%.2f, phi=%.1f°, psi=%.1f°, v=[%.2f,%.2f,%.2f]\n',...
                state.rx, state.rz, state.phi, state.psi, state.vx, state.vy, state.vz);
    end
end

function accepted = accept_solution(current_cost, new_cost, T)
    delta = new_cost - current_cost;
    if delta < 0
        accepted = true;
    else
        accepted = rand() < exp(-delta/T);
    end
end

function [dlambda, dmu1, dmu2, dphi, dpsi] = apply_constraints(dl, dm1, dm2, dp, dpsi, state)
    % 发动机开度变化率约束 (防止突变)
    dlambda = max(-0.05, min(0.05, dl));  % 限制开度变化率在±5%/s
    
    % 摆角变化率约束
    max_ang_rate = 2; % °/s
    dmu1 = max(-max_ang_rate, min(max_ang_rate, dm1));
    dmu2 = max(-max_ang_rate, min(max_ang_rate, dm2));
    
    % 滚转/偏航变化率约束
    max_rot_rate = 5; % °/s
    dphi = max(-max_rot_rate, min(max_rot_rate, dp));
    dpsi = max(-max_rot_rate, min(max_rot_rate, dpsi));
    
    % 结合当前状态的二次约束
    % 示例：如果当前phi接近下限，限制减少方向的变化
    if state.phi < 65 && dphi < 0
        dphi = max(dphi, 0.5); % 在接近下限时禁止继续减小
    end
end

%% 保存数据
function save_trajectory_data(state, controls, dt)
    % 生成时间序列
    time = (0:length(state)-1)*dt;
    
    % 转换状态数据为结构体数组
    data_struct = struct(...
        'time', num2cell(time'),...
        'rx', num2cell([state.rx]'),...
        'ry', num2cell([state.ry]'),...
        'rz', num2cell([state.rz]'),...
        'vx', num2cell([state.vx]'),...
        'vy', num2cell([state.vy]'),...
        'vz', num2cell([state.vz]'),...
        'm', num2cell([state.m]'),...
        'lambda', num2cell([state.lambda]'),...
        'mu1', num2cell([state.mu1]'),...
        'mu2', num2cell([state.mu2]'),...
        'phi', num2cell([state.phi]'),...
        'psi', num2cell([state.psi]'),...
        'dlambda', num2cell(controls(1,:)'),...
        'dmu1', num2cell(controls(2,:)'),...
        'dmu2', num2cell(controls(3,:)'),...
        'dphi', num2cell(controls(4,:)'),...
        'dpsi', num2cell(controls(5,:)')...
    );

    % 保存为MAT文件
    save('trajectory_data.mat', 'data_struct', '-v7.3');
    
    % 另存为CSV文件（可选）
%     csv_data = [
%         time', 
%         [state.rx]', [state.ry]', [state.rz]',
%         [state.vx]', [state.vy]', [state.vz]',
%         [state.m]', [state.lambda]',
%         [state.mu1]', [state.mu2]',
%         [state.phi]', [state.psi]',
%         controls'
%     ];
%     
%     header = {
%         'Time(s)', 'Xpos(m)', 'Ypos(m)', 'Zpos(m)',...
%         'Vx(m/s)', 'Vy(m/s)', 'Vz(m/s)', 'Mass(kg)','lambda',...
%         'mu1(rad)', 'mu2(rad)', 'phi(rad)', 'psi(rad)',...
%         'dlambda', 'dmu1', 'dmu2', 'dphi', 'dpsi'
%     };
%     
%     fid = fopen('trajectory_data.csv', 'w');
%     fprintf(fid, '%s,', header{1:end-1});
%     fprintf(fid, '%s\n', header{end});
%     fclose(fid);
%     
%     dlmwrite('trajectory_data.csv', csv_data, '-append', 'precision', '%.6f');
%     
    disp('轨迹数据已保存为:');
    disp(['   MAT文件: ' fullfile(pwd, 'trajectory_data.mat')]);
%     disp(['   CSV文件: ' fullfile(pwd, 'trajectory_data.csv')]);
end

%% 可视化函数
function visualize_results(state, controls, dt)
    time = 0:dt:(length(state)-1)*dt;
    
%     figure('Position', [100,100,1200,800])
    
    % 三维轨迹
%     subplot(2,3,1)
    figure;
    plot3([state.rx], [state.rz], [state.ry], 'k-o')
    xlabel('x'); ylabel('z'); zlabel('y');
    grid on
    print -dpng 三维.png
    
    time = time +414;
    % 速度分量
    figure;
    plot(time, [state.rx], 'r-o', time, [state.ry], 'g-o', time, [state.rz], 'b-o')
    xlabel('时间(s)');
    legend('x轴坐标值(m)','y轴坐标值(m)','z轴坐标值(m)');
    grid on
    print -dpng 位置.png
    
    % 速度分量
    figure;
    plot(time, [state.vx], 'r-o', time, [state.vy], 'g-o', time, [state.vz], 'b-o')
    xlabel('时间(s)');
    legend('x轴速度值(m/s)','y轴速度值(m/s)','z轴速度值(m/s)');
    grid on
    print -dpng 速度.png

    % 质量变化
    figure;
    plot(time, [state.m], 'r-o', 'LineWidth',2)
    xlabel('时间(s)'); ylabel('质量(kg)');
    grid on
    print -dpng 质量.png

    % 发动机开度
    figure;
    plot(time, [state.lambda], 'r-o')
    xlabel('时间(s)');ylabel('发动机开度');
    grid on
    print -dpng 发动机开度.png

    % 发动机摆角
    figure;
    plot(time, [state.mu1], 'r-o',...
         time, [state.mu2], 'b-o')
    xlabel('时间(s)');
    legend('发动机摆角1(°)','发动机摆角2(°)');
    grid on
    print -dpng 发动机摆角.png

    % 姿态角
    figure;
    plot(time, [state.phi], 'r-o',...
         time, [state.psi], 'b-o')
    xlabel('时间(s)');
    legend('俯仰角(°)','偏航角(°)');
    grid on
    print -dpng 姿态角.png

    % 发动机开度变化
    figure;
    plot(time, controls(1,:), 'r-o');
    xlabel('时间(s)');ylabel('发动机开度变化率');
    grid on
    print -dpng 发动机开度变化.png

    % 发动机摆角变化
    figure;
    plot(time, controls(2,:), 'r-o', time, controls(3,:), 'b-o')
    xlabel('时间(s)');
    legend('发动机摆角1变化率(°/s)','发动机摆角2变化率(°/s)');
    grid on
    print -dpng 发动机摆角变化.png

    % 姿态角变化
    figure;
    plot(time, controls(4,:), 'r-o',time, controls(5,:), 'b-o');
    xlabel('时间(s)');
    legend('俯仰角变化率(°/s)','偏航角变化率(°/s)');
    grid on
    print -dpng 姿态角变化.png

    % 终端状态显示
    fprintf('终端质量: %.2f kg\n', state(end).m);
    fprintf('着陆位置: %.2f m, %.2f m, %.2f m\n', norm([state(end).rx, state(end).ry, state(end).rz]));
    fprintf('着陆速度: %.2f m/s, %.2f m/s, %.2f m/s\n', norm([state(end).vx, state(end).vy, state(end).vz]));
end