
% a horrid accreted mess to draw a channelflow schematic diagram
% little rhyme or reason
clf;
hold;

% define size of box
Lx = 2;
Ly = 1;
Lz = 1;

% define shrinking factors for perspective of viewpoint
py = 0.4;
pyLz = py*Lz;
px = 0.6;
pxLz = px*Lz;

grey = [.9 .9 .9];

axis([ -0.5 3 -1 2.5])

% draw opaque rectangles for walls
patch([0 Lx Lx+pxLz pxLz], [0 0 pyLz pyLz], grey)
patch([0 Lx Lx+pxLz pxLz], Ly+[0 0 pyLz pyLz], grey)

% draw edges of box, include hidden edge with dots
plot([0 0], [0 Ly], 'k');
plot([Lx Lx], [0 Ly], 'k');
plot([Lx+pxLz Lx+pxLz], [pyLz pyLz+Ly], 'k');
plot([pxLz pxLz], [pyLz Ly], 'k');
plot([pxLz pxLz], [Ly pyLz+Ly], 'k:');

% draw mean velocity profile
N=50;
y = 0:Ly/N:Ly;
yp = 2*y - Ly;
U = (1.1-yp.^4 - 0.1*yp.^2)/2;
%U = (y-y.^2);
%plot(Lz/2 + [Lx/2 Lx/2], pyLz/2+ [0 Ly], 'k')
%plot(Lz/2 + Lx/2+ 2*(y-y.^2), pyLz/2+y, 'k')
a = Lx/2;
b = 0;
plot([a a], [0 Ly], 'k')
plot(a+ U, y, 'k')

% draw arrows within mean velocity profile
N = 10;
y = Ly/N:Ly/N:Ly-Ly/N;
yp = 2*y - Ly;
U = (1.1-yp.^4- 0.1*yp.^2)/2;
%U = (y-y.^2);
d4 = 0.08;
d2 = 0.04;
d1 = 0.02;
for n = 1 :N-1
  plot([a+U(n) a], [y(n) y(n)], 'k');
  plot([a+U(n) a+U(n)-d2], [y(n) y(n)-d1], 'k');
  plot([a+U(n) a+U(n)-d2], [y(n) y(n)+d1], 'k');
end
text(a+max(U), b+0.75*Ly, 'U(y)')

% draw coordinate axes
zx = 0;   % x-position of origin
zy = -0.3;  % y-position of origin
Lxp = 0.2*Lx;
Lyp = 0.6*Ly;
Lzp = 0.6*Lz;
Ly2 = Ly/2;

% the positioning of labels is fudged manually, with d1,d2,d4
% adjustments and visual review of results

plot([zx zx], [zy zy+Lyp], 'k');
plot([zx zx-d1], [zy+Lyp zy+Lyp-d2], 'k');
plot([zx zx+d1], [zy+Lyp zy+Lyp-d2], 'k');
text(zx+d2+d1, zy+Lyp, 'y');

plot([zx zx+Lxp], [zy zy], 'k');
plot([zx+Lxp zx+Lxp-d2], [zy zy-d1], 'k');
plot([zx+Lxp zx+Lxp-d2], [zy zy+d1], 'k');
text(zx+Lxp, zy-0.1, 'x');

plot([zx zx+px*Lzp], [zy zy+py*Lzp], 'k');
plot([zx+px*Lzp zx+px*Lzp-d1 ], [zy+py*Lzp zy+py*Lzp-d2], 'k');
plot([zx+px*Lzp zx+px*Lzp-d2 ], [zy+py*Lzp zy+py*Lzp], 'k');
text(zx+px*Lzp+d2, zy+py*Lzp-d1, 'z');


text(zx-d4, zy-d2, '0');

% 
zx = 0;
zy = Ly/2;

text(zx-0.1, Ly, 'b');
text(zx-0.1, 0, 'a');
text(zx+Lx/2+px*Lz*2/3, Ly+py*Lz+0.07, 'L_x');
text(zx+Lx+px*Lz/2+0.05, py*Lz/2-0.1, 'L_z');
axis off
print -deps2 schematic.eps