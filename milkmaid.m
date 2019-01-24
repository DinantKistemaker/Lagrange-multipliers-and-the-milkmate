% Cool example of using the optimization technique "Lagrange Multipliers".
% Here we solve the milkmaid problem to illustrate this technique.
% This is part of Introduction to Optimal Control, in the bachelor course
% "Biophysics of Locomotion", taught at the VU University Amsterdam. This
% technique will later be used to optimize dynamical musculoskeletal
% models...
%
%
% The milkmaid is at position mm and needs to go to the river before heading back
% to the farm at location fa. She has lots of milk to carry and is in a
% hurry, and thus she wants to go to the farm with the shortest route possible
% while still going to the river. The river streams according to:
% x2 = (1/6*x1.^3 + .5*x1.^2 - x1 +5). Which path should she take?
% This is a constraint optimization problem and is solved here using Lagrange 
% multipliers together with Newton-Rapson's method of root finding.
% L = f(x) +lambda*c(x); with L the Lagrangian, f(x) the cost function and 
% c(x) the constraint function. Extremum of the Lagrangian is found when
% the partial derivative of L w.r.t both x and lambda are zero. These
% necessary conditions are transformed into a KKT system using first order
% Taylor approximations. All derrivatives are calculated using the symbolic
% toolbox.
%
% Dinant Kistemaker 21-01-2016. dinant.kistemaker@gmail.com

clear all, close all, clc
%% plotting the problem
x1=-10:.1:10;
x2=1/6*x1.^3 + .5*x1.^2 - x1 +5;
mm=[-2;4];
fa=[0;0];

X0 = [2 2]; %Initial guess of optimal point

hold on
plot(x1,x2,'b','Linewidth',6)
plot(mm(1),mm(2),'ro','MarkerFaceColor','r'),plot(fa(1),fa(2),'ko','MarkerFaceColor','k')
plot(X0(1),X0(2),'bo','MarkerFaceColor','b')
axis([-6 6 -2 10]), axis equal
l=legend('river','milkmaid','farm','initial (dumb) guess');
l.Location='SouthEast';
l.AutoUpdate = 'off';
xlabel('x'),ylabel('y')
set(gca,'FontName','Times','FontSize',18)
% print -dtiff -r600 milmaid_problem

x2=-10:.1:10;
for i=1:length(x1)
    for j=1:length(x2)
       f(i,j)=sqrt((x1(i)-mm(1)).^2+(x2(j)-mm(2)).^2) + sqrt((fa(1)-x1(i)).^2+(fa(2)-x2(j)).^2);
    end
end
contour(x1,x2,f',20)
% print -dtiff -r600 milmaid_contour
disp('Done plotting problem...')

%% Doing the actual Optimization
clearallbut X0
syms x1 x2 lambda
disp('Milkmaid position')
mm=[-2;4]
pause
disp('Farm position')
fa=[0;0]
pause
disp('River position = constraint')
c= x2 -(1/6*x1.^3 + .5*x1.^2 - x1 +5);
pretty(c),pause
disp('distance from the milkmaid to a point [x,y] and then to the farm = cost')
f=sqrt((x1-mm(1))^2+(x2-mm(2))^2) + sqrt((fa(1)-x1)^2+(fa(2)-x2)^2);
pretty(f),pause
disp('Jacobian of constraint')
J=[diff(c,x1) diff(c,x2)];
pretty(J),pause
disp('gradient of cost')
g=[diff(f,x1);diff(f,x2)];
pretty(g),pause

disp('Hessian of constraint')
Hc=[diff(J(1),x1) diff(J(1),x2);diff(J(2),x1) diff(J(2),x2)];
pretty(Hc),pause
disp('Hessian of cost')
Hf=[diff(g,x1) diff(g,x2)];
pretty(Hf);,pause
disp('Hessian of Lagrangian')
Hl=Hf-lambda*Hc;
pretty(Hl),pause

disp('KKT system')
A = [Hl -J';J 0];
pretty(A),pause
b = [-g+J'*lambda;-c];
pretty(b),pause

disp('Initial gues')
x1=X0(1)
x2=X0(2)
lambda=1
X0=[x1 x2 lambda]';

disp('Hit any key to start Newton Iterations')
pause

max_round=200;
c_tol=1e-4;
x_tol=1e-4;
for round=1:max_round
    plot(x1,x2,'go','MarkerFaceColor','g')
    U=eval(A);
    s=eval(b);
    dX = U\s
    x1    =x1    +dX(1);
    x2    =x2    +dX(2);
    lambda=lambda+dX(3);
    nr_iteration_rounds=round-1
    if (eval(abs(c))<c_tol) & (abs(dX)<x_tol)
        disp('Necessary conditions reached within set tollerance')
        plot(x1,x2,'yo','MarkerFaceColor','y')
        break
    end
    pause
    
end
