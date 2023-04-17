clear all

global N; N = 2;   %number of generalised coordinates
syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
global l1 l2 m1 m2; l1 = 1; l2 = 1; m1 = 1; m2 = 1;   %vars for double pendulum

"Our Lagrangian is:"
L = makeLagrangian(1,0)

"Lagrange equations in generalised coordinates and velocities:"
Eqs = lagrangeEquations(L)

"Lagrange equation system as derivatives of generalised coordinates:"
Sys = lagrangeEquationsSolvable(L)

"Can we use dsolve to solve this system?"
Sol = dsolve(Sys == zeros(N,1))

"ODE45 solution to this system of equations:"
M = 10;  %number of swings
plotSystem(Sys,[0:0.01:M*pi],[pi 0 pi -0.0001]);

function L = makeLagrangian(g,c)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global l1 l2 m1 m2; %double pendulum params

    T = 0.5*m1*l1*l1*u1*u1 + 0.5*m2*(l1*l1*u1*u1 + l2*l2*u2*u2 + 2*l1*l2*u1*u2*cos(q1-q2)); %kinetic energy
    V = -(m1 + m2)*g*l1*cos(q1) - m2*g*l2*cos(q2); %potential energy
    L = exp(c*t)*(T - V);   %simulated damping with the exponential term

end

function Eqs = lagrangeEquations(L)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global l1 l2 m1 m2; %double pendulum params

    Eqs = [];

    for i = 1:N

        dLdqi = diff(L, str2sym("q" + string(i) + "(t)"));  %differentiate by qi
        dLdui = diff(L, str2sym("u" + string(i) + "(t)"));  %differentiate by ui
        dt = diff(dLdui,t); %differentiate by time

        for j = 1:N

            dt = subs(dt, diff(str2sym("q" + string(j) + "(t)"),t), str2sym("u" + string(j) + "(t)"));  %replace differentiated qi's with ui's

        end

        Eqs = [Eqs; dt - dLdqi]; %lagrange equation for i'th generalised coordinate

    end

end

function Sys = lagrangeEquationsSolvable(L)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global l1 l2 m1 m2; %double pendulum params

    Sys = [];

    for i = 1:N

        dLdqi = diff(L, str2sym("q" + string(i) + "(t)"));  %differentiate by qi
        dLdui = diff(L, str2sym("u" + string(i) + "(t)"));  %differentiate by ui
        dt = diff(dLdui,t); %differentiate by time
        dX = dt - dLdqi;

        for j = 1:N

            dX = subs(dX, diff(str2sym("u" + string(j) + "(t)"),t), diff(str2sym("q" + string(j) + "(t)"),t,2)); %replace differentiated ui's with 2 * differentiaed qi's
            dX = subs(dX, str2sym("u" + string(j) + "(t)"), diff(str2sym("q" + string(j) + "(t)"),t)); %replace ui's with differentiated qi's

        end

        Sys = [Sys; dX]; %lagrange equation for i'th generalised coordinate

    end

end

function plotSystem(Sys, tspan, q0)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global l1 l2 m1 m2; %double pendulum params
%     hold on;

    eq = odeToVectorField(Sys == zeros(N,1));
    odefun = matlabFunction(eq, 'Vars', {'t','Y'});

    [t,Q] = ode45(odefun, tspan, q0);

    %plot(t,Q(:,[1:2:2*N]));

    animateSystem(t,Q);

end

function animateSystem(t,Q)

    global l1 l2 m1 m2; %double pendulum params
    figure;

    P = [l1*sin(Q(:,1)) + l2*sin(Q(:,3)) -l1*cos(Q(:,1)) - l2*cos(Q(:,3))]; %pendulum end
    R = [l1*sin(Q(:,1)) -l1*cos(Q(:,1))]; %pendulum middle

    xlim([-(l1+l2) l1+l2]);
    ylim([-(l1+l2) l1+l2]);
    
    h = animatedline('Marker','o','Color','red');
    hold on;
    
    P1 = plot([0 0]);
    P2 = plot(P(1,:) - R(1,:));
    P3 = plot(R(1,:));
    P4 = plot(R(1,:));

    %P(3,1)  %gets 3rd item in first angle values

    for j = 1:length(P(:,1))
       clearpoints(h);
       addpoints(h, P(j,1), P(j,2));
       delete(P1);
       delete(P2);
       delete(P3);
       delete(P4);
       P1 = [P1 plot([0 R(j,1)],[0 R(j,2)], 'Color', 'black')]; %draws first arm
       P2 = [P1 plot([R(j,1) P(j,1)],[R(j,2) P(j,2)], 'Color', 'black')];   %draws second arm
       P3 = [P2 plot(R([1:j],1), R([1:j],2), 'Color', 'green')]; %draws curve of first arm bob
       P4 = [P2 plot(P([1:j],1), P([1:j],2), 'Color', 'blue')]; %draws curve of second arm bob
       drawnow;

    end

end