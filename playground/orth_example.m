function orth_example()
clear all
close all

x = 4 + 1i;
y = 1 + 2*1i;

%x = x/abs(x);
%y = y/abs(y);

compass(real(x),imag(x),'b');
hold
compass(real(x),imag(x'),'b--');

compass(real(y),imag(y),'r');

x_coord = [real(x) imag(x)];
y_coord = [real(y) imag(y)];



[x_orth_coord , proj_y] = get_ortho(x_coord,y_coord)
[y_orth_coord ,proj_x ] = get_ortho(y_coord,x_coord)

compass(x_orth_coord(1),x_orth_coord(2),'g');

compass(y_orth_coord(1),y_orth_coord(2),'k');


compass(proj_x(1),proj_x(2),'m');

compass(proj_y(1),proj_y(2),'c');

x_ortho2 = get_ortho2(x,y);

compass(real(x_ortho2),imag(x_ortho2),'y');

function  [ortho, proj ] = get_ortho(x,y)
    
    proj  = ((x*y')/(y*y'))*y;
    ortho = x - proj;
    
function ortho = get_ortho2(x,y)
    
    ortho = x*conj(y)/abs(y);