function [alt, vel] = VelocityUKF(z, dt)
%
%
persistent Q R
persistent x P
persistent init

if isempty(init)
    Q = [ 0  0;
          0  0.0018 ];

    R = 10;

    x = [ 1.5  0 ]';
    P = 5 * eye(2);

    init = 1;
end

[Xi, W] = SigmaPoints(x, P, 0);

%fXi = zeros(n, 2 * n + 1);
fXi = zeros(2, 5);
for k = 1 : 5
    fXi(:, k) = fx(Xi(:, k), dt);
end

[xp, Pp] = UT(fXi, W, Q);


%hXi = zeros(m, 2 * n + 1);
hXi = zeros(1, 5);
%for k = 1 : 2 * n + 1
for k = 1 : 5
    hXi(:, k) = hx(fXi(:, k));
end

[zp, Pz] = UT(hXi, W, R);

%Pxz = zeros(n, m);
Pxz = zeros(2, 1);
%for k = 1 : 2 * n + 1
for k = 1 : 5
    Pxz = Pxz + W(k) * (fXi(:, k) - xp) * (hXi(:, k) - zp)';
end

K = Pxz * inv(Pz);

x = xp + K * (z - zp);
P = Pp - K * Pz * K';

alt = x(1);
vel = x(2);


%-----------------------------
function xp = fx(x, dt)
%
%
A = [ 1  dt;
      0  1 ];
xp = A * x;
%-----------------------------
function yp = hx(x)
%
%
yp = [ 1  0 ] * x;