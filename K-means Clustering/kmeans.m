%kmeans.m
%Richard Chen

k2 = fopen('k2.dat', 'r');
k3 = fopen('k3.dat', 'r');
k4 = fopen('k4.dat', 'r');
k5 = fopen('k5.dat', 'r');
k6 = fopen('k6.dat', 'r');

k2_data = textscan(k2, '%f %f %f %d');
k3_data = textscan(k3, '%f %f %f %d');
k4_data = textscan(k4, '%f %f %f %d');
k5_data = textscan(k5, '%f %f %f %d');
k6_data = textscan(k6, '%f %f %f %d');

k2_col_x = k2_data{1};
k2_col_y = k2_data{2};
k2_col_z = k2_data{3};
k2_col_c = k2_data{4};

k3_col_x = k3_data{1};
k3_col_y = k3_data{2};
k3_col_z = k3_data{3};
k3_col_c = k3_data{4};

k4_col_x = k4_data{1};
k4_col_y = k4_data{2};
k4_col_z = k4_data{3};
k4_col_c = k4_data{4};

k5_col_x = k5_data{1};
k5_col_y = k5_data{2};
k5_col_z = k5_data{3};
k5_col_c = k5_data{4};

k6_col_x = k6_data{1};
k6_col_y = k6_data{2};
k6_col_z = k6_data{3};
k6_col_c = k6_data{4};

figure(1);
hold on;
scatter3(k2_col_x, k2_col_y, k2_col_z, 10, k2_col_c, 'filled')
hold off;

figure(2);
hold on;
scatter3(k3_col_x, k3_col_y, k3_col_z, 10, k3_col_c, 'filled')
hold off;

figure(3);
hold on;
scatter3(k4_col_x, k4_col_y, k4_col_z, 10, k4_col_c, 'filled')
hold off;

figure(4);
hold on;
scatter3(k5_col_x, k5_col_y, k5_col_z, 10, k5_col_c, 'filled')
hold off;

figure(5);
hold on;
scatter3(k6_col_x, k6_col_y, k6_col_z, 10, k6_col_c, 'filled')
hold off;