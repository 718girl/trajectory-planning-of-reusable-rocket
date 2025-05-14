function [state_history, control_history] = SA_rocket_trajectory()
    tic;
    %% ������ʼ��
    N = 50;            % ʱ��ڵ���
    dt = 0.8;            % ʱ�䲽��(s)
    T0 = 1000;         % ��ʼ�¶�
    alpha = 0.95;      % ����ϵ��
    max_iter = 50;     % ÿ�¶ȵ�������
    
    %% ���Ʊ����߽� (d��/dt, d��1/dt, d��2/dt, d��/dt, d��/dt)
    control_lim = [-0.3,  0.3;       % d��/dt 
                  -5, 5;      % d��1/dt [��/s]
                  -5, 5;      % d��2/dt [��/s]
                  -5, 5;    % d��/dt [��/s]
                  -5, 5]';  % d��/dt [��/s]
    
    %% ģ���˻��ʼ��
    current_controls = initialize_controls(N, control_lim);
    [current_cost, current_state] = evaluate_trajectory(current_controls, dt);
    
    %% �˻���ѭ��
    T = T0;
    while T > 1e-3
        for iter = 1:max_iter
            % ���������
            new_controls = generate_neighbor(current_controls, control_lim, T/T0);
            
            % �����켣
            [new_cost, new_state] = evaluate_trajectory(new_controls, dt);
            
            % Metropolis׼��
            if accept_solution(current_cost, new_cost, T)
                current_controls = new_controls;
                current_cost = new_cost;
                current_state = new_state;
            end
        end
        T = T * alpha;
    end
    
    %% ������
    total_time = toc;
    fprintf('������ʱ��Ϊ%.5f\n', total_time);
    [cost, state_terminal, state_history, control_history] = evaluate_trajectory(current_controls, dt);
    fprintf('���ۺ������%.5f\n', cost);
    save_trajectory_data(state_history, control_history, dt)
    visualize_results(state_history, control_history, dt);
end

%% �ؼ��Ӻ���
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
    
    % ʩ�ӱ߽�Լ��
    for i = 1:size(lim,2)
        new_controls(i,:) = max(lim(1,i), min(lim(2,i), new_controls(i,:)));
    end
end

function [cost, state, state_history, control_history] = evaluate_trajectory(controls, dt)
    % ״̬������ʼ��
    state.rx = -600;        % Xλ��(m)
    state.ry = 6000;        % Y�߶�(m)
    state.rz = -200;        % Zλ��(m)
    state.vx = -80;         % X�ٶ�(m/s)
    state.vy = -300;        % Y�ٶ�(m/s)
    state.vz = 30;          % Z�ٶ�(m/s)
    state.m = 49000;        % ����(kg)
    state.lambda = 0.35;    % ���������ȣ���ʼ��Ϊ���ޣ�
    state.mu1 = 5;          % �������ڽ�(��)
    state.mu2 = 5;          % �������ڽ�(��)
    state.phi = 80;         % ������(��)
    state.psi = -3;         % ƫ����(��)
    
    n_steps = size(controls,2);
    state_history = repmat(state, [1 n_steps]);
    control_history = zeros(5, n_steps);
    penalty = 0; % ����Լ���ͷ�

    for k = 1:n_steps
        % ================== ����Լ������ ==================
        % Լ��1: ����������lambda �� [0.35,1]��ӲԼ����
        state.lambda = max(0.35, min(1, state.lambda));
        
        % Լ��2: �������ڽ�mu1,mu2 �� [-5,5]�㣨ӲԼ����
        state.mu1 = max(-5, min(5, state.mu1));
        state.mu2 = max(-5, min(5, state.mu2));
        
        % Լ��3: phi �� [60,90]�㣨ӲԼ����
        state.phi = max(60, min(90, state.phi));
        
        % Լ��4: psi �� [-5,5]�㣨ӲԼ����
        state.psi = max(-5, min(5, state.psi));
        
        % Լ��5: ����m �� 28000kg��ӲԼ����
        if state.m < 28000
            cost = inf;
            return; % ������ֹ����
        end
        
        % Լ��6: �߶ȵ����ݼ�����Լ����
        if k > 1 && state.ry > state_history(k-1).ry
            penalty = penalty + 1e5 * (state.ry - state_history(k-1).ry);
        end
        
        % ================== ���������� ==================
        raw_dlambda = controls(1,k);
        raw_dmu1 = controls(2,k);
        raw_dmu2 = controls(3,k);
        raw_dphi = controls(4,k);
        raw_dpsi = controls(5,k);
        
        % ================== ����Լ������ ==================
        % Լ�������Ŀ�����
        [dlambda, dmu1, dmu2, dphi, dpsi] = apply_constraints(...
            raw_dlambda, raw_dmu1, raw_dmu2, raw_dphi, raw_dpsi, state);
        
        % ��¼ʵ��Ӧ�õĿ�����
        control_history(:,k) = [dlambda; dmu1; dmu2; dphi; dpsi]; % �洢Լ����Ŀ�����
        
        % ���¿������״̬����ʩ��ӲԼ����
        state.lambda = state.lambda + dlambda*dt;
        state.mu1 = state.mu1 + dmu1*dt;
        state.mu2 = state.mu2 + dmu2*dt;
        state.phi = state.phi + dphi*dt;
        state.psi = state.psi + dpsi*dt;
        
        % ================== ����ѧ���� ==================
        Tmax = 912000 - 67000 * exp(-state.ry/7110);
        T = state.lambda * Tmax;
        
        % �����������㣨ʹ�ýǶ�Լ�����ֵ��
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
        
        % ״̬����
        state.rx = state.rx + state.vx*dt;
        state.ry = state.ry + state.vy*dt;
        state.rz = state.rz + state.vz*dt;
        state.vx = state.vx + (Fx + Dx)/state.m * dt;
        state.vy = state.vy + (Fy + Dy)/state.m * dt - 9.807*dt;
        state.vz = state.vz + (Fz + Dz)/state.m * dt;
        state.m = state.m - T/(isp*9.807)*dt;
        
        state_history(k) = state;
    end
    
    % ================== �ն�Լ������ ==================
    terminal_penalty = 0;
    
    % �߶�Լ�� ry=18.6��0.1m
    if abs(state.ry - 18.6) > 0.1
        terminal_penalty = terminal_penalty + 1e8*(state.ry - 18.6)^2; %8 
    end
    
    % λ��Լ�� rx,rz �� [-25,25]m
    if abs(state.rx) > 25 || abs(state.rz) >25
        terminal_penalty = terminal_penalty + 1e8*(max(abs(state.rx),abs(state.rz))-25)^2; %7 ��֤�ն�λ��8
    end
    
    % ��̬Լ�� phi=90��1��, psi=0��0.5��
    if abs(state.phi - 90) > 1
        terminal_penalty = terminal_penalty + 1e10*(state.phi - 90)^2;
    end
    if abs(state.psi) > 0.5
        terminal_penalty = terminal_penalty + 1e10*(state.psi)^2;
    end
    
    % �ٶ�Լ�� vx,vz��[-3,3], vy��[-1,1] m/s
    if abs(state.vx) >3 || abs(state.vz)>3
        terminal_penalty = terminal_penalty + 1e9*(max(abs(state.vx),abs(state.vz))-3)^2;  %10 ��֤�ٶ�11
    end
    if abs(state.vy) >1
        terminal_penalty = terminal_penalty + 1e9*(abs(state.vy)-1)^2;
    end
    
    % �ۺϴ��ۼ���
    cost = -state.m + penalty + terminal_penalty;
    
    % ���ӻ���飨����ʱ���ã�
    if rand < 0.01 % 1%���ʳ������
        fprintf('�ն�״̬���: rx=%.2f, rz=%.2f, phi=%.1f��, psi=%.1f��, v=[%.2f,%.2f,%.2f]\n',...
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
    % ���������ȱ仯��Լ�� (��ֹͻ��)
    dlambda = max(-0.05, min(0.05, dl));  % ���ƿ��ȱ仯���ڡ�5%/s
    
    % �ڽǱ仯��Լ��
    max_ang_rate = 2; % ��/s
    dmu1 = max(-max_ang_rate, min(max_ang_rate, dm1));
    dmu2 = max(-max_ang_rate, min(max_ang_rate, dm2));
    
    % ��ת/ƫ���仯��Լ��
    max_rot_rate = 5; % ��/s
    dphi = max(-max_rot_rate, min(max_rot_rate, dp));
    dpsi = max(-max_rot_rate, min(max_rot_rate, dpsi));
    
    % ��ϵ�ǰ״̬�Ķ���Լ��
    % ʾ���������ǰphi�ӽ����ޣ����Ƽ��ٷ���ı仯
    if state.phi < 65 && dphi < 0
        dphi = max(dphi, 0.5); % �ڽӽ�����ʱ��ֹ������С
    end
end

%% ��������
function save_trajectory_data(state, controls, dt)
    % ����ʱ������
    time = (0:length(state)-1)*dt;
    
    % ת��״̬����Ϊ�ṹ������
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

    % ����ΪMAT�ļ�
    save('trajectory_data.mat', 'data_struct', '-v7.3');
    
    % ���ΪCSV�ļ�����ѡ��
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
    disp('�켣�����ѱ���Ϊ:');
    disp(['   MAT�ļ�: ' fullfile(pwd, 'trajectory_data.mat')]);
%     disp(['   CSV�ļ�: ' fullfile(pwd, 'trajectory_data.csv')]);
end

%% ���ӻ�����
function visualize_results(state, controls, dt)
    time = 0:dt:(length(state)-1)*dt;
    
%     figure('Position', [100,100,1200,800])
    
    % ��ά�켣
%     subplot(2,3,1)
    figure;
    plot3([state.rx], [state.rz], [state.ry], 'k-o')
    xlabel('x'); ylabel('z'); zlabel('y');
    grid on
    print -dpng ��ά.png
    
    time = time +414;
    % �ٶȷ���
    figure;
    plot(time, [state.rx], 'r-o', time, [state.ry], 'g-o', time, [state.rz], 'b-o')
    xlabel('ʱ��(s)');
    legend('x������ֵ(m)','y������ֵ(m)','z������ֵ(m)');
    grid on
    print -dpng λ��.png
    
    % �ٶȷ���
    figure;
    plot(time, [state.vx], 'r-o', time, [state.vy], 'g-o', time, [state.vz], 'b-o')
    xlabel('ʱ��(s)');
    legend('x���ٶ�ֵ(m/s)','y���ٶ�ֵ(m/s)','z���ٶ�ֵ(m/s)');
    grid on
    print -dpng �ٶ�.png

    % �����仯
    figure;
    plot(time, [state.m], 'r-o', 'LineWidth',2)
    xlabel('ʱ��(s)'); ylabel('����(kg)');
    grid on
    print -dpng ����.png

    % ����������
    figure;
    plot(time, [state.lambda], 'r-o')
    xlabel('ʱ��(s)');ylabel('����������');
    grid on
    print -dpng ����������.png

    % �������ڽ�
    figure;
    plot(time, [state.mu1], 'r-o',...
         time, [state.mu2], 'b-o')
    xlabel('ʱ��(s)');
    legend('�������ڽ�1(��)','�������ڽ�2(��)');
    grid on
    print -dpng �������ڽ�.png

    % ��̬��
    figure;
    plot(time, [state.phi], 'r-o',...
         time, [state.psi], 'b-o')
    xlabel('ʱ��(s)');
    legend('������(��)','ƫ����(��)');
    grid on
    print -dpng ��̬��.png

    % ���������ȱ仯
    figure;
    plot(time, controls(1,:), 'r-o');
    xlabel('ʱ��(s)');ylabel('���������ȱ仯��');
    grid on
    print -dpng ���������ȱ仯.png

    % �������ڽǱ仯
    figure;
    plot(time, controls(2,:), 'r-o', time, controls(3,:), 'b-o')
    xlabel('ʱ��(s)');
    legend('�������ڽ�1�仯��(��/s)','�������ڽ�2�仯��(��/s)');
    grid on
    print -dpng �������ڽǱ仯.png

    % ��̬�Ǳ仯
    figure;
    plot(time, controls(4,:), 'r-o',time, controls(5,:), 'b-o');
    xlabel('ʱ��(s)');
    legend('�����Ǳ仯��(��/s)','ƫ���Ǳ仯��(��/s)');
    grid on
    print -dpng ��̬�Ǳ仯.png

    % �ն�״̬��ʾ
    fprintf('�ն�����: %.2f kg\n', state(end).m);
    fprintf('��½λ��: %.2f m, %.2f m, %.2f m\n', norm([state(end).rx, state(end).ry, state(end).rz]));
    fprintf('��½�ٶ�: %.2f m/s, %.2f m/s, %.2f m/s\n', norm([state(end).vx, state(end).vy, state(end).vz]));
end