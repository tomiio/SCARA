function create_cylinder(x,y,z,r,d,color)
%r = 0.75; d = 3;x = 2;y=2;z=2;color = 'red';
[X0,Y0,Z0] = cylinder(r,100);
surf(X0+x,Y0+y,Z0*d+z,'facecolor',color,'LineStyle','none');
hold on;
fill3(X0(1,:) + x, Y0(1,:) + y, Z0(1,:) + z, color);
fill3(X0(2,:) + x, Y0(2,:) + y, Z0(2,:)*d + z, color);