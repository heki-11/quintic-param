clear all
close all
%% Create parameter functions in l,xu,yu,xuu,yuu
load('x_but')
load('y_but')
totlen = length(x_but);
flag = 0;
ind = 1;
inddif = 4;
ind_accum= [];
xknot= [];
yknot= [];
while ind < length(x_but)
    ind = ind + inddif;
    if ind > length(x_but)
        break
    end
    xknot = [xknot;x_but(ind)];
    yknot = [yknot;y_but(ind)];
    ind_accum = [ind_accum;ind];
    if ind == 1
        ind = 2;
    end
end
save('ind_accum.mat','ind_accum')
xknot = [xknot;x_but(end)];
yknot = [yknot;y_but(end)];
no_point = length(xknot);
save('no_point.mat','no_point')
figure;plot(xknot,yknot,'bo')
axis equal;drawnow;
save('xknot.mat','xknot')
save('yknot.mat','yknot')
len = length(xknot)
syms l [len-1 1]
syms xu [len 1] 
syms yu [len 1] 
syms xuu [len 1] 
syms yuu [len 1]
syms u [len-1 1]
syms xuk
syms yuk
syms J
syms ObjFunc
for i = 1:len-1
    Ax(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^5*((6*(xknot(i+1)-xknot(i)))-3*(xu(i+1)+xu(i))*l(i)+1/2*(xuu(i+1)-xuu(i))*l(i)^2);
    Bx(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^4*((15*(xknot(i)-xknot(i+1)))+(7*xu(i+1)+8*xu(i))*l(i)+(3/2*xuu(i)-xuu(i+1))*l(i)^2);
    Cx(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^3*((10*(xknot(i+1)-xknot(i)))-(4*xu(i+1)+6*xu(i))*l(i)-(3/2*xuu(i)-1/2*xuu(i+1))*l(i)^2);
    Dx(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/2*xuu(i);
    Ex(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        xu(i);
    
    Ay(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^5*((6*(yknot(i+1)-yknot(i)))-3*(yu(i+1)+yu(i))*l(i)+1/2*(yuu(i+1)-yuu(i))*l(i)^2);
    By(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^4*((15*(yknot(i)-yknot(i+1)))+(7*yu(i+1)+8*yu(i))*l(i)+(3/2*yuu(i)-yuu(i+1))*l(i)^2);
    Cy(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/l(i)^3*((10*(yknot(i+1)-yknot(i)))-(4*yu(i+1)+6*yu(i))*l(i)-(3/2*yuu(i)-1/2*yuu(i+1))*l(i)^2);
    Dy(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        1/2*yuu(i);
    Ey(l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        yu(i);
    
    xuk(u(i),l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        5*Ax*u(i)^4 + 4*Bx*u(i)^3 + 3*Cx*u(i)^2 + 2*Dx*u(i) + Ex;
    yuk(u(i),l(i),xu(i+1),yu(i+1),xuu(i+1),yuu(i+1),xu(i),yu(i),xuu(i),yuu(i)) = ...
        5*Ay*u(i)^4 + 4*By*u(i)^3 + 3*Cy*u(i)^2 + 2*Dy*u(i) + Ey;
    J = (xuk^2+yuk^2-1)^2;
    ObjFunc(i) = int(J,u(i),0,l(i));
    i
    
end
save('ObjFunc.mat','ObjFunc')
clear
%% Solve for ds/du-optimal Quintic Spline
load('ObjFunc')
load('xknot')
load('yknot')
load('x_but')
load('y_but')
totlen = length(x_but);
load('ind_accum')
len = length(xknot);
fun_scalar = sum(ObjFunc);
syms l [len-1 1]
syms xu [len 1] 
syms yu [len 1] 
syms xuu [len 1] 
syms yuu [len 1]
x = [l;xu;xuu;yu;yuu];
fun = matlabFunction(fun_scalar,'vars',{x});
arclen = zeros(1,length(xknot)-1);
load('no_point')
for i = 1:length(ind_accum)-1
    if i < length(ind_accum)-1
        x_incr = [x_but(ind_accum(i):ind_accum(i+1))]';
        y_incr = [y_but(ind_accum(i):ind_accum(i+1))]';
    else
        x_incr = [x_but(ind_accum(i):end)]';
        y_incr = [y_but(ind_accum(i):end)]';
    end
    for j = 2:length(x_incr)
        arclen(i) = arclen(i) + sqrt((x_incr(j)-x_incr(j-1))^2+(y_incr(j)-y_incr(j-1))^2);
    end
end
arclen = [arclen(1:end-1) arclen(2) arclen(1)];
save('arclen.mat','arclen')
arclen_accum = [cumsum(arclen)];
xu0 = gradient(xknot,arclen_accum);
yu0 = gradient(yknot,arclen_accum);
xuu0 = gradient(xu0,arclen_accum);
yuu0 = gradient(yu0,arclen_accum);
x0 = [[arclen(1:end-1)]';xu0;xuu0;yu0;yuu0];
options = optimoptions('fmincon','display','iter','maxfunctionevaluations',1e+8,...
    'algorithm','sqp',...
    'maxiterations',100,'plotfcn','optimplotfval');
Aineq = sparse(-[[eye(len-1) zeros(len-1,4*len)]]);
bineq = sparse([zeros(len-1,1)]);

[sol,fval,exitflag,output] = fmincon(fun,x0,Aineq,bineq,[],[],[],[],[],options);
save('sol4.mat','sol')

figure;
plot(1:len-1,x0(1:len-1),'bo-');
hold on;
plot(1:len-1,sol(1:len-1),'rx-');
legend('Initial guess of l','Optimal solution of l')

figure;
subplot(2,1,1)
plot(1:len,x0(len:2*len-1),'bo-');
hold on;
plot(1:len,sol(len:2*len-1),'rx-');
legend('Initial guess of xu','Optimal solution of xu')

subplot(2,1,2)
plot(1:len,x0(len:2*len-1),'bo-');
hold on;
plot(1:len,sol(len:2*len-1),'rx-');
legend('Initial guess of yu','Optimal solution of yu')

%% Retrieve parameters
l = sol(1:len-1);
xu = sol(len:2*len-1);
xuu = sol(2*len:3*len-1);
yu = sol(3*len:4*len-1);
yuu = sol(4*len:5*len-1);
load('arclen')
arclen_accum = cumsum(arclen);
l2 = arclen';
xu2 = gradient(xknot,arclen_accum);
yu2 = gradient(yknot,arclen_accum);
xuu2 = gradient(xu2,arclen_accum);
yuu2 = gradient(yu2,arclen_accum);
sol2 = [arclen';xu2;xuu2;yu2;yuu2];
save('sol2.mat','sol2')

for i = 1:len-1
    Ax(i) = 1/l(i)^5*((6*(xknot(i+1)-xknot(i)))-3*(xu(i+1)+xu(i))*l(i)+1/2*(xuu(i+1)-xuu(i))*l(i)^2);
    Bx(i) = 1/l(i)^4*((15*(xknot(i)-xknot(i+1)))+(7*xu(i+1)+8*xu(i))*l(i)+(3/2*xuu(i)-xuu(i+1))*l(i)^2);
    Cx(i) = 1/l(i)^3*((10*(xknot(i+1)-xknot(i)))-(4*xu(i+1)+6*xu(i))*l(i)-(3/2*xuu(i)-1/2*xuu(i+1))*l(i)^2);
    Dx(i) = 1/2*xuu(i);
    Ex(i) = xu(i);
    Fx(i) = xknot(i);
    
    Ay(i) = 1/l(i)^5*((6*(yknot(i+1)-yknot(i)))-3*(yu(i+1)+yu(i))*l(i)+1/2*(yuu(i+1)-yuu(i))*l(i)^2);
    By(i) = 1/l(i)^4*((15*(yknot(i)-yknot(i+1)))+(7*yu(i+1)+8*yu(i))*l(i)+(3/2*yuu(i)-yuu(i+1))*l(i)^2);
    Cy(i) = 1/l(i)^3*((10*(yknot(i+1)-yknot(i)))-(4*yu(i+1)+6*yu(i))*l(i)-(3/2*yuu(i)-1/2*yuu(i+1))*l(i)^2);
    Dy(i) = 1/2*yuu(i);
    Ey(i) = yu(i);
    Fy(i) = yknot(i);
    
    Ax2(i) = 1/l2(i)^5*((6*(xknot(i+1)-xknot(i)))-3*(xu2(i+1)+xu2(i))*l2(i)+1/2*(xuu2(i+1)-xuu2(i))*l2(i)^2);
    Bx2(i) = 1/l2(i)^4*((15*(xknot(i)-xknot(i+1)))+(7*xu2(i+1)+8*xu2(i))*l2(i)+(3/2*xuu2(i)-xuu2(i+1))*l2(i)^2);
    Cx2(i) = 1/l2(i)^3*((10*(xknot(i+1)-xknot(i)))-(4*xu2(i+1)+6*xu2(i))*l2(i)-(3/2*xuu2(i)-1/2*xuu2(i+1))*l2(i)^2);
    Dx2(i) = 1/2*xuu2(i);
    Ex2(i) = xu2(i);
    Fx2(i) = xknot(i);
    
    Ay2(i) = 1/l2(i)^5*((6*(yknot(i+1)-yknot(i)))-3*(yu2(i+1)+yu2(i))*l2(i)+1/2*(yuu2(i+1)-yuu2(i))*l2(i)^2);
    By2(i) = 1/l2(i)^4*((15*(yknot(i)-yknot(i+1)))+(7*yu2(i+1)+8*yu2(i))*l2(i)+(3/2*yuu2(i)-yuu2(i+1))*l2(i)^2);
    Cy2(i) = 1/l2(i)^3*((10*(yknot(i+1)-yknot(i)))-(4*yu2(i+1)+6*yu2(i))*l2(i)-(3/2*yuu2(i)-1/2*yuu2(i+1))*l2(i)^2);
    Dy2(i) = 1/2*yuu2(i);
    Ey2(i) = yu2(i);
    Fy2(i) = yknot(i);

end

QuinticCoef = [Ax;Bx;Cx;Dx;Ex;Fx;Ay;By;Cy;Dy;Ey;Fy];
QuinticCoef2 = [Ax2;Bx2;Cx2;Dx2;Ex2;Fx2;Ay2;By2;Cy2;Dy2;Ey2;Fy2];
save('QuinticCoef4.mat','QuinticCoef')
save('QuinticCoef2.mat','QuinticCoef2')

for i = 1:len-1
    param = [0:0.001*l(i):l(i)];
    xcurve = Ax(i)*param.^5 + Bx(i)*param.^4 + Cx(i)*param.^3 + Dx(i)*param.^2 + Ex(i)*param + Fx(i);
    ycurve = Ay(i)*param.^5 + By(i)*param.^4 + Cy(i)*param.^3 + Dy(i)*param.^2 + Ey(i)*param + Fy(i);
    curve(:,2*i-1:2*i) = [xcurve' ycurve'];
end

xcurve_tot = [];
ycurve_tot = [];

for j = 1:len-1
    xcurve_tot = [xcurve_tot;curve(:,2*j-1)];
    ycurve_tot = [ycurve_tot;curve(:,2*j)];
end

load('x_but')
load('y_but')

figure;
plot(xcurve_tot,ycurve_tot,'b','linewidth',2);
hold on
plot(x_but,y_but,'r--','linewidth',2)
legend('ds/du-optimal curve','original curve')
