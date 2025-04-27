function [alt, vel] = VelocityUKFWithAccZ(z, dt, a)

a = (a - 1) * 9.80665;

persistent Q R x P init
if isempty(init)
    Q = [0.001, 0;
         0,    0.00001];
    R = 10;

    x = [1.5; 0];

    P = [2.0, 0;
         0, 1.0];

    init = true;
end

[Xi, W] = SigmaPoints(x, P, 0);
L = numel(W);

fXi = zeros(2, L);
for k = 1:L
    fXi(:,k) = fx(Xi(:,k), dt, a);
end
[xp, Pp] = UT(fXi, W, Q);

if ~isempty(z)
    hXi = zeros(1, L);
    for k = 1:L
        hXi(k) = hx(fXi(:,k));
    end
    [zp, Pz] = UT(hXi, W, R);
    Pxz = zeros(2,1);
    for k = 1:L
        Pxz = Pxz + W(k) * (fXi(:,k) - xp) * (hXi(k) - zp);
    end
    K = Pxz * inv(Pz);
    
    x = xp + K * (z - zp);
    P = Pp - K * Pz * K';
else
    x = xp;
    P = Pp;
end

alt = x(1);
vel = x(2);
end

%-----------------------------
function xp = fx(x, dt, a)
    A = [1, dt;
         0, 1];
    xp = A * x + [0.5 * a * dt^2; a * dt];
end

%-----------------------------
function yp = hx(x)
    yp = x(1);
end