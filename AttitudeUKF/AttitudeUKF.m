function q = AttitudeUKF(z_q, dt, omega_meas)

    persistent Q R x P init
    if isempty(init)
        Q = eye(4) * 1e-6;

        R = eye(4) * 1e-3;

        x = [1 0 0 0]';

        P = eye(4) * 0.1;
        
        init = true;
    end
    
    %Make column vectors
    if ~isempty(z_q)
        z_q = z_q(:);
    end
    omega = omega_meas(:);
    
    %state vector dimension and sigma points
    n = 4;
    L = 2*n + 1;

    [X, W] = SigmaPoints(x, P, 0);   % X:4×L, W:1×L
    
    %Prediction
    Xp = zeros(n, L);
    for k = 1:L
        %predict sigma points and integrates quaternion
        qk = X(:,k);
        q_pred = fx_quat(qk, dt, omega);
        %maintaining unit size by normalization
        Xp(:,k) = q_pred / norm(q_pred);
    end
    [x_pred, P_pred] = UT(Xp, W, Q);
    x_pred = x_pred(:);

    if ~isempty(z_q)
        Hq = Xp;
        [z_pred, Pz] = UT(Hq, W, R);
        z_pred = z_pred(:);

        Pxz = zeros(n, n);
        for k = 1:L
            Pxz = Pxz + W(k) * (Xp(:,k) - x_pred) * (Hq(:,k) - z_pred)';
        end
        %Kalman gain
        K = Pxz / Pz;
        
        %column vector
        innov = z_q - z_pred;
        innov = innov(:);
        
        %update state and covariance
        x = x_pred + K * innov;
        x = x / norm(x);
        P = P_pred - K * Pz * K';
    else
        x = x_pred;
        P = P_pred;
    end

    q = x;
end

function q_out = fx_quat(q, dt, omega)
    omega = omega(:);
    Omega = [ 0      -omega';
              omega  -skew(omega) ];
    q_out = q + 0.5 * Omega * q * dt;
end

function S = skew(v)
    v = v(:);
    S = [  0    -v(3)  v(2);
          v(3)   0    -v(1);
         -v(2)  v(1)   0   ];
end