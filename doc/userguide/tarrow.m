% arrow(x0, y0, x3, y3, h)
% draw arrow with tip at x0,y0, head at x1,y1, headsize=h
function t = tarrow(x0, y0, x3, y3, h, lbl);

phi = atan2(y3-y0, x3-x0);
c = cos(phi);
s = sin(phi);
h2 = h/2;
h4 = h/4;

% arrowhead base pt
x1 = x3 - h*c;
y1 = y3 - h*s;

% arrowhead corner
x2 = x1 + h4*s;
y2 = y1 - h4*c;

% arrowhead corner
x4 = x1 - h4*s;
y4 = y1 + h4*c;

% 
x5 = x0 + h4*s;
y5 = y0 - h4*c;
x6 = x0 - h4*s;
y6 = y0 + h4*c;

% label position
xl = x3 + h*c;
yl = y3 + h*c;

x = [x1 x2 x3 x4 x1 x0 x5 x6] ;
y = [y1 y2 y3 y4 y1 y0 y5 y6];
plot(x,y,'k');

text(xl,yl,lbl, 'FontSize',16);

