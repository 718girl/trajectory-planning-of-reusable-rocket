function gpops_rocket_reentry_full()
% ȫ�ֲ���
m0 = 49000;                         % ��ʼ����(kg)
m_stru = 28000;                     % �ṹ����(kg)
lambda_min2 = 0.35;                 % ��������������
lambda_max = 1.0;                   % ��������������
lambda_rmax = 0.3;                  % ���������ȱ仯������޶�
mu_max = 5;                         % �������ڽ����Ƕ�(��)
mu_rmax = 5;                        % �������ڽǵ����仯��(��/s)
w_max3 = 5;                         % ���������ת������3(��/s)

% �߽硢����״̬Լ��
bounds.phase.initialtime.lower = 414;
bounds.phase.initialtime.upper = 414;
bounds.phase.initialstate.lower = [-600, 6000, -200, -80, -300, 30, m0, 0.8, 5, 5, 80, -3];
bounds.phase.initialstate.upper = [-600, 6000, -200, -80, -300, 30, m0, 0.8, 5, 5, 80, -3];

bounds.phase.state.lower =  [-1000, 18.6, -200, -200, -500, -100, m_stru, lambda_min2,     0,      0,   60, -5];
bounds.phase.state.upper =  [100,   6000,  100,  200,   10,  100,     m0, lambda_max, mu_max, mu_max,   100, 10];
bounds.phase.control.lower =  [-lambda_rmax, -mu_rmax, -mu_rmax, -w_max3, -w_max3];
bounds.phase.control.upper =  [ lambda_rmax,  mu_rmax,  mu_rmax,  w_max3,  w_max3];

% ����ʽ(3.25)�ն�Լ��
bounds.phase.finaltime.lower = 440;      
bounds.phase.finaltime.upper = 480;
bounds.phase.finalstate.lower = [-25, 18.6, -25, -3, -1, -3, m_stru, lambda_min2,      0,      0,  89.9, -0.1]; 
bounds.phase.finalstate.upper = [ 25, 18.6,  25,  3,  1,  3,     m0,  lambda_max, mu_max, mu_max,  90.1,  0.1]; % ��ֱ��̬90��0.1��

% ��ʼ�²�ֵ

data = load('guess.mat'); % �������ɵĲ²�����
% data = load('initinal_guess.mat'); % �������ɵĲ²�����
phase1 = data.phase6;

guess.phase.time = phase1.time;
guess.phase.state = phase1.states;
guess.phase.control = phase1.controls;

%  GPOPS-II������
setup.name = 'RLV_Full_Trajectory';
setup.functions.continuous = @rocket_dynamics;
setup.functions.endpoint = @endpoint;
setup.bounds = bounds;
setup.guess = guess; 
setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.mesh.method = 'hp1';    % hp-PattersonRao
setup.mesh.tolerance = 1e-6;
setup.mesh.maxiterations = 5;
setup.mesh.phase.colpoints = 50;  % * ones(1, 10)
setup.mesh.phase.fraction = 0.1;  %  * ones(1, 10)

% Ŀ�꺯��
function output = endpoint(input)
    output.objective = -input.phase.finalstate(7);
end

totalTic = tic;
output = gpops2(setup);
totalTime = toc(totalTic);
disp(totalTime);

% ����ѧ����
function phase_dynamics = rocket_dynamics(input)
    phase = input.phase; 
    t = phase.time;
    x = phase.state;
    u = phase.control;

    rx = x(:,1); ry = x(:,2); rz = x(:,3);
    vx = x(:,4); vy = x(:,5); vz = x(:,6);
    m = x(:,7); lambda = x(:,8); 
    mu1 = x(:,9); mu2 = x(:,10);
    phi = x(:,11); psi = x(:,12);

    dlambda = u(:,1); 
    dmu1 = u(:,2); dmu2 = u(:,3);
    omega_z = u(:,4); omega_y = u(:,5);
    
    % ����������
    Tmax = 912000 - 67000 * exp(-ry/7110);
    T = lambda.*Tmax;

    % ��������ϵ��������
    Tx_body = T .* cosd(mu1) .* cosd(mu2);
    Ty_body = T .* cosd(mu1) .* sind(mu2);
    Tz_body = -T .* sind(mu1);

    % ת��Ϊ��½����ϵ
%     Fx = Tx_body .* cosd(phi) .* cosd(psi) - Ty_body .* sind(phi) .* cosd(psi) + Tz_body .* sind(psi);
%     Fy = Tx_body .* sind(phi) + Ty_body .* cosd(phi);
%     Fz = -Tx_body .* cosd(phi) .* sind(psi) + Ty_body .* sind(phi) .* sind(psi) - Tz_body .* cosd(psi);
    Fx = Tx_body .* cosd(phi) .* cosd(psi) - Ty_body .* sind(phi) + Tz_body .*cosd(phi).* sind(psi);
    Fy = Tx_body .* sind(phi) .* cosd(psi) + Ty_body .* cosd(phi) + Tz_body .*sind(phi).* sind(psi);
    Fz = -Tx_body .* sind(psi) + Tz_body .* cosd(psi);

    % ��������(��½����ϵ)
    rho = 1.225 * exp(-ry/7110);
    V = sqrt(vx.^2 + vy.^2 + vz.^2);
    Dx = -0.5 .* rho .* 10.52 * 1.5 .* V .* vx;
    Dy = -0.5 .* rho .* 10.52 * 1.5 .* V .* vy;
    Dz = -0.5 .* rho .* 10.52 * 1.5 .* V .* vz;
    
    isp = 312 - 29 * exp(-ry/7110);
    
    % ����ѧ���̣���½����ϵ��
    drx = vx;
    dry = vy;
    drz = vz;
    dvx = (Fx + Dx) ./ m; 
    dvy = (Fy + Dy) ./ m - 9.807;
    dvz = (Fz + Dz) ./ m;  
    dm = -T ./ (isp * 9.807);  
    dphi = omega_z;
    dpsi = omega_y;

    phase_dynamics.dynamics = [drx, dry, drz, dvx, dvy, dvz, dm, dlambda, dmu1, dmu2, dphi, dpsi];
end

%% ��ͼ

solution = output.result.solution;
phase1 = solution.phase;
save( 'reference_trajectory.mat', 'phase1');

t6 = solution.phase.time;   % ��ȷ����λ

r6 = solution.phase.state(:,1:3);
rx = r6(:,1);
ry = r6(:,2);
rz = r6(:,3);
v6 = solution.phase.state(:,4:6);
m6 = solution.phase.state(:,7);
lambda6 = solution.phase.state(:,8);
mu16 = solution.phase.state(:,9);
mu26 = solution.phase.state(:,10);
phi6 = solution.phase.state(:,11);
psi6 = solution.phase.state(:,12);

dlambda6 = solution.phase.control(:,1);
dmu16 = solution.phase.control(:,2);
dmu26 = solution.phase.control(:,3);
dphi6 = solution.phase.control(:,4);
dpsi6 = solution.phase.control(:,5);

% ��ά�켣����
figure;
plot3(rx, rz, ry, 'k-', 'LineWidth', 2);
grid on;
hold on;
plot3(rx(1), rz(1), ry(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot3(rx(end), rz(end), ry(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
set(gca, 'FontSize', 14);
xlabel('x');
ylabel('z');
zlabel('y');
print -dpng ��άͼ.png

% x
figure;
plot(t6, r6(:, 1), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('x������ֵ(m)');
print -dpng x.png

% y
figure;
plot(t6, r6(:, 2), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('y������ֵ(m)');
print -dpng y.png

% z
figure;
plot(t6, r6(:, 3), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('z������ֵ(m)');
print -dpng z.png

% vx
figure;
plot(t6, v6(:, 1), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('x���ٶ�ֵ(m/s)');
print -dpng dx.png

% vy
figure;
plot(t6, v6(:, 2), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('y���ٶ�ֵ(m/s)');
print -dpng dy.png

% vz
figure;
plot(t6, v6(:, 3), 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('z���ٶ�ֵ(m/s)');
print -dpng dz.png

% m
figure;
plot(t6, m6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('����(kg)');
print -dpng m.png

% lambda
figure;
plot(t6, lambda6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('����������');
print -dpng lambda.png

% mu1
figure;
plot(t6, mu16, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('�������ڽ�1(��)');
print -dpng mu1.png

% mu2
figure;
plot(t6, mu26, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('�������ڽ�2(��)');
print -dpng mu2.png

% phi
figure;
plot(t6, phi6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('������(��)');
print -dpng phi.png

% psi
figure;
plot(t6, psi6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('ƫ����(��)');
print -dpng psi.png

% dlambda
figure;
plot(t6, dlambda6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('���������ȱ仯��');
print -dpng dlambda.png

% dmu1
figure;
plot(t6, dmu16, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('�������ڽ�1�仯��');
print -dpng dmu1.png

% dmu2
figure;
plot(t6, dmu26, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('�������ڽ�2�仯��');
print -dpng dmu2.png

% dphi
figure;
plot(t6, dphi6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('�����Ǳ仯��');
print -dpng dphi.png

% dpsi
figure;
plot(t6, dpsi6, 'r-', 'LineWidth', 1.5);
grid on;
set(gca, 'FontSize', 14);
xlabel('ʱ��(s)');
ylabel('ƫ���Ǳ仯��');
print -dpng dpsi.png
end