function [alt, vel] = VelocityUKF(z, dt, a)
% VelocityUKF: UKF for altitude and vertical velocity
%   z:    (optional) measured altitude [m]
%   dt:   sampling interval [s]
%   a:    measured vertical acceleration [m/s^2]
%   alt, vel: estimated altitude and velocity

persistent Q R x P init
if isempty(init)
    % Process noise covariance for [z; v]
    Q = [0.01, 0;
         0,    0.1];
    % Measurement noise covariance (altitude)
    R = 10;
    % Initial state [z0; v0]
    x = [1.5; 0];
    % Initial covariance
    P = 5 * eye(2);
    init = true;
end

% Generate sigma points (n=2 -> 2n+1 points)
[Xi, W] = SigmaPoints(x, P, 0);
L = numel(W);

% Prediction step: propagate each sigma-point via fx
fXi = zeros(2, L);
for k = 1:L
    fXi(:,k) = fx(Xi(:,k), dt, a);
end
[xp, Pp] = UT(fXi, W, Q);

% Update step if measurement available
if ~isempty(z)
    % Measurement sigma-points via hx
    hXi = zeros(1, L);
    for k = 1:L
        hXi(k) = hx(fXi(:,k));
    end
    [zp, Pz] = UT(hXi, W, R);
    % Cross-covariance Pxz (2x1)
    Pxz = zeros(2,1);
    for k = 1:L
        Pxz = Pxz + W(k) * (fXi(:,k) - xp) * (hXi(k) - zp);
    end
    % Kalman gain
    K = Pxz / Pz;
    % State and covariance update
    x = xp + K * (z - zp);
    P = Pp - K * Pz * K';
else
    % No measurement: carry forward prediction
    x = xp;
    P = Pp;
end

% Extract estimates
alt = x(1);
vel = x(2);
end

%-----------------------------
function xp = fx(x, dt, a)
% Process model: discrete kinematics with acceleration input
    A = [1, dt;
         0, 1];
    xp = A * x + [0.5 * a * dt^2; a * dt];
end

%-----------------------------
function yp = hx(x)
% Measurement model: extract altitude component
    yp = x(1);
end

%-----------------------------
%% Test Script
% To run these tests, clear the function state: clear VelocityUKF
fprintf('--- VelocityUKF Test Cases ---\n');

dt = 1; a_val = 0;
clear VelocityUKF;  % reset persistent
[alt1, vel1] = VelocityUKF([], dt, a_val);
fprintf('Test1 (no input, a=0): alt=%.3f, vel=%.3f\n', alt1, vel1);

% Constant acceleration test
dt = 1; a_val = 1;
fprintf('\nConstant accel (a=1 m/s^2)\n');
clear VelocityUKF;
for k = 1:5
    [alt_k, vel_k] = VelocityUKF([], dt, a_val);
    fprintf(' Step %d: alt=%.3f, vel=%.3f\n', k, alt_k, vel_k);
end

% Measurement update test (true height = 10)
dt = 1; a_val = 0;
fprintf('\nMeasurement update (z=10)\n');
clear VelocityUKF;
for k = 1:5
    [alt_k, vel_k] = VelocityUKF(10, dt, a_val);
    fprintf(' Update %d: alt=%.3f, vel=%.3f\n', k, alt_k, vel_k);
end