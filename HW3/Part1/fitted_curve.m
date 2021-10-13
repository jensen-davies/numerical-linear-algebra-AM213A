clf
curve = importdata('curve_Chol.dat');
atkinson = importdata('atkinson.dat');
x_atkinson = atkinson.data(:,1);
y_atkinson = atkinson.data(:,2);

hold on
scatter(x_atkinson, y_atkinson);
plot(x_atkinson,curve);
title('Cholesky Degree 3 Fitted Curve of Atkinson Data');
xlabel('Atkinson x');
ylabel('Atkinson y');
legend('Atkinson Data','Fitted Curve');

