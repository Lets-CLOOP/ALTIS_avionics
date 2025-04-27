clear all

dt = 0.002;
t  = 0 : dt : 65.245;

sampleCounts = length(t);

estimated = zeros(sampleCounts, 2);
predicted = zeros(sampleCounts, 1);

for k = 1 : sampleCounts
  z = GetAlt();      
  [alt, vel] = VelocityUKF(z, dt);

  estimated(k, :) = [alt, vel];
  predicted(k)    = z;
end


figure
hold on
title('Filtered altitude')
plot(t, estimated(:, 1), 'b.')
xlabel('time(s)')
ylabel('altitude(m)')

figure
hold on
title('Calculated velocity')
plot(t, estimated(:, 2), 'b.')
xlabel('time(s)')
ylabel('velocity(m/s)')
