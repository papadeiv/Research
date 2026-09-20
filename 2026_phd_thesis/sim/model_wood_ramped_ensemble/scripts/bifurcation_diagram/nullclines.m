H = -0.38;
z = @(x) BoxModel_2DH_IVP (0, x, H, 'FamousB_2xCO2');

figure; hold on;

contour(X, Y, F1, [0 0], 'r', 'LineWidth', 2);  % f1 = 0
contour(X, Y, F2, [0 0], 'b', 'LineWidth', 2);  % f2 = 0

legend('f_1 = 0', 'f_2 = 0');
xlabel('x'); ylabel('y');
title('Nullclines');
grid on;