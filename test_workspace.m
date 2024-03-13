x = 0; y = 0; z = 0;
[X0,Y0,Z0] = cylinder(3,1000);
opacity = 0.2;

p = surf(X0(:,1:347) + x, Y0(:,1:347) + y, Z0(:,1:347)*1.6 + z,"EdgeColor","none");
p.FaceAlpha = opacity;
hold on;
p = surf(X0(:,654:1001) + x, Y0(:,654:1001) + y, Z0(:,654:1001)*1.6 + z,"EdgeColor","none");
p.FaceAlpha = opacity;

t = linspace(-125*pi/180,125*pi/180);
a = 1.25; b = 1.75;
x0=0; y0=0;              % circles center
c = sqrt(a^2 + b^2 -2*a*b*cos(pi-145*pi/180));
rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
rout = 3;   % radii sizes
% edge = linspace(
% edge2 = 
p =fill3([x0+rout*cos(t),x0+flip(rin*cos(t))],[y0+rout*sin(t),y0+flip(rin*sin(t))],(1.6)*ones(1,2*length(rout*cos(t))),'r','linestyle','none');
p.FaceAlpha = opacity;
p =fill3([x0+rout*cos(t),x0+flip(rin*cos(t))],[y0+rout*sin(t),y0+flip(rin*sin(t))],(0)*ones(1,2*length(rout*cos(t))),'r','linestyle','none');
p.FaceAlpha = opacity;

c = sqrt(a^2 + b^2 -2*a*b*cos(pi-145*pi/180));
[X1,Y1,Z1] = cylinder(c,1000);
p = surf(X1 + x, Y1 + y, Z1*1.6 + z,"EdgeColor","none");
p.FaceAlpha = opacity;

y = sin(55*pi/180)*1.25; x = -cos(55*pi/180)*1.25; z = 0;
h = sin(55*pi/180)*1.25;
alpha = acos(h/1.75) + (35*pi/180);
beta = pi - alpha;
[X0,Y0,Z0] = cylinder(1.75,1000);
j = 601;
i = 345;
p = surf(X0(:,i:j) + x, Y0(:,i:j) + y, Z0(:,i:j)*1.6 + z,"EdgeColor","none");
p.FaceAlpha = opacity;

t = linspace(125*pi/180,beta+125*pi/180);
t1 = linspace(125*pi/180,pi);
x0=-sqrt(1.25^2 - h^2); y0=h;              % circles center
rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
rout = 1.75;   % radii sizes
p = fill3([x0+rout*cos(t),0+flip(rin*cos(t1))],[y0+rout*sin(t),0+flip(rin*sin(t1))],(1.6)*ones(1,2*length(rout*cos(t))),'red','linestyle','none');
p.FaceAlpha = opacity;
p = fill3([x0+rout*cos(t),0+flip(rin*cos(t1))],[y0+rout*sin(t),0+flip(rin*sin(t1))],(0)*ones(1,2*length(rout*cos(t))),'red','linestyle','none');
p.FaceAlpha = opacity;

y = -y;
[X0,Y0,Z0] = cylinder(1.75,1000);
j = 670;
i = 400;
p = surf(X0(:,i:j) + x, Y0(:,i:j) + y, Z0(:,i:j)*1.6 + z,"EdgeColor","none");
p.FaceAlpha = opacity;

t = linspace(alpha+55*pi/180,alpha+beta+55*pi/180);
t1 = linspace(pi,alpha+beta+55*pi/180);
a = 1.25; b = 1.75;
x0=-sqrt(1.25^2 - h^2); y0= -h;              % circles center
rin = c;%sqrt(1.25^2 + 1.75^2 - 2*1.25*1.75*cos(pi-2)); 
rout = 1.75;   % radii sizes
tx = t;
p = fill3([x0+rout*cos(tx),0+flip(rin*cos(t1))],[y0+rout*sin(tx),0+flip(rin*sin(t1))],(1.6)*ones(1,2*length(rout*cos(tx))),'red','linestyle','none');
p.FaceAlpha = opacity;
p = fill3([x0+rout*cos(tx),0+flip(rin*cos(t1))],[y0+rout*sin(tx),0+flip(rin*sin(t1))],(0)*ones(1,2*length(rout*cos(tx))),'red','linestyle','none');
p.FaceAlpha = opacity;

