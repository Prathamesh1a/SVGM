clc
clear
% Import acceleration vs time data
data=readmatrix("Time-History.xlsx");%Reads input time history at reference site
acceleration = data(:,2);
time = data(:,1);

% Calculate squared acceleration
acceleration_squared = acceleration.^2;

% Integrate squared acceleration vs time data
integral = cumtrapz(time, acceleration_squared);

% Calculate cumulative count of zero-level up-crossings
zero_crossings = find(acceleration(1:end-1).*acceleration(2:end) < 0);
up_crossings = zero_crossings(acceleration(zero_crossings+1) >= 0);
cumulative_up_crossings = cumsum(ones(size(up_crossings)));

% Calculate cumulative count of positive minima
[minima, min_locs] = findpeaks(-acceleration);
positive_minima = -minima(minima < 0);
cumulative_positive_minima = cumsum(ismember(-minima, positive_minima));

% Create time vectors for minima and maxima
min_time = time(min_locs);

% Plot the results
figure;
subplot(4,1,1);
plot(time, acceleration);
xlabel('Time');
ylabel('Acceleration');
title('Acceleration vs Time');

subplot(4,1,2);
plot(time, integral);
xlabel('Time');
ylabel('Integral of Squared Acceleration');
title('Integral of Squared Acceleration vs Time');

subplot(4,1,3);
plot(time(up_crossings), cumulative_up_crossings);
xlabel('Time');
ylabel('Cumulative Count');
title('Cumulative Count of Zero-Level Up-Crossings');

subplot(4,1,4);
plot(min_time(ismember(-minima, positive_minima)), cumulative_positive_minima(ismember(-minima, positive_minima)), 'g');
xlabel('Time');
ylabel('Cumulative Count');
title('Cumulative Count of Positive Minima');
