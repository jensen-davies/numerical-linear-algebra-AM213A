[x,y] = meshgrid(-5:0.5:5);

a = 3.90335731207141490E-002;
b = 0.36338249407750045;
c = 7.80671462414282979E-002;

z = 1/c - (a/c)*x - (b/c)*y;

clf; hold on;
surf(x,y,z, 'DisplayName','Plane','FaceColor',[0 0.4470 0.7410], 'EdgeColor','none')

plot3(1,2,3,'k', 'marker','o','DisplayName','A','MarkerFaceColor','r')
plot3(-3,2,5,'k','marker','o','DisplayName','B', 'MarkerFaceColor','y')
plot3(pi, exp(1), -sqrt(2), 'k','marker','o','DisplayName','C', 'MarkerFaceColor','m')

legend

