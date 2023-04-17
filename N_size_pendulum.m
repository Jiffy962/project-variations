clear vars

%number of pendulum bobs. we need this just about everywhere
global N; N = 3;

%syms vars for generalised coordinates and velocities
%qi's are anticlockwise angles of each bob from the vertical, and ui's are anticlockwise angular velocities of each bob
syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);

%main vars
global li mi;  %lengths of arms, and masses of bobs per arm

%li = [1 1 1 1]; mi = [1 1 1 1]; %specify custom values for li, mi (will have to append for larger values of N)
li = ones(1,N); mi = ones(1,N);    %each bob has equal mass/length

%q0 is a vector indicating initial conditions for the generalised coordinates and velocities
%odd pos are for initial generalised coordinates qi_0, even values are for initial generalised velocities ui_0

%q0 = [0.5*pi 0 0.5*pi 0 0.5*pi 0 0.5*pi 0 0.5*pi 0] %specify custom values for q0 (will have to append for larger values of N)
q0 = createInitialConditions(); %generate initial conditions where the pendulum is balanced, with a small nudge for the top bob

tspan = [0:0.005:40*pi];  %timespan to run ode45 for

%end of program/vars prep
%main program

%generates a Lagrangian for N pendulums - if second argument > 0, an artificial damping is applied
%note that changing g=9.8 only changes how fast the simulation runs, as the potential is proportional to g
"Our Lagrangian is:"
L = makeLagrangian(9.8,0)

%generates Euler-Lagrange equations by differentiating Lagrangian by generalised coordinates and velocities
"Lagrange equations in generalised coordinates and velocities:"
Eqs = lagrangeEquations(L)

%generates Euler-Lagrange equations by differentiating Lagrangian by generalised coordinates and velocities
%Sys is then usable for functions like dsolve, ode45, etc
"Lagrange equation system as derivatives of generalised coordinates:"
Sys = lagrangeEquationsSolvable(L)

%if system can be solved directly, we have a printed solution
%this rarely ever works if we do not apply taylor series approximations for sine/cosine
"Can we use dsolve to solve this system?"
Sol = dsolve(Sys == zeros(N,1))

%solves system by converting to first-order vector field and applying ode45
"Running ODE45 solution to this system of equations:"
[t,Q] = solveSystem(Sys,tspan,q0);

%plots solved system as an animated swinging pendulum
animateSystem(t,Q);

%end of main program
%member functions

function q0 = createInitialConditions()

    global N;

    A = [pi 0];
    q0 = repmat(A,1,N);
    q0(2*N) = -0.01;    %small nudge

end

function L = makeLagrangian(g,c)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global li mi; %double pendulum params

    %prealloc
    T = 0; V = 0; xi = 0; yi = 0; dxi = 0; dyi = 0;
    
    for i = 1:N
            
        qi = str2sym("q" + string(i) + "(t)");
        ui = str2sym("u" + string(i) + "(t)");

        %x/y/dx/dy components
        xi = xi + li(i)*sin(qi);
        yi = yi - li(i)*cos(qi);
        dxi = dxi + li(i) * ui * cos(qi);
        dyi = dyi + li(i) * ui * sin(qi);

        T = T + 0.5*mi(i)*(dxi*dxi + dyi*dyi);   %kinetic energy
        V = V + mi(i)*g*yi;  %potential energy

    end

    L = exp(c*t)*(T - V);   %simulated damping with the exponential term

end

function Eqs = lagrangeEquations(L)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities

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

function [t,Q] = solveSystem(Sys, tspan, q0)

    global N; syms t; syms q(t) [1 N]; syms u(t) [1 N]; u(t) = diff(q(t),1,t);    %generalised coordinates and velocities
    global li mi; %double pendulum params

    eq = odeToVectorField(Sys == zeros(N,1));
    odefun = matlabFunction(eq, 'Vars', {'t','Y'});

    [t,Q] = ode45(odefun, tspan, q0);

end

function animateSystem(t,Q)

    global N li; %double pendulum params
    figure;

    plotsize = sum(li);
    steps = length(Q(:,1));
    L = 100;
    
    xlim([-plotsize plotsize]);
    ylim([-plotsize plotsize]);
    pbaspect([1 1 1]);
    
    h = animatedline('Marker','o','Color','black');
    hold on;
    grid on;

    %prealloc
    P = [];
    xi = 0; yi = 0;
    bobs(:,:,1) = [zeros(steps,1) zeros(steps,1)];

    %calculate points of each bob on the pendulum; i represents the i'th bob from the centre
    %to reference: bobs(timestep,x=1/y=2,i)
    for i = 1:N

        xi = xi + li(i)*sin(Q(:,2*i - 1));
        yi = yi - li(i)*cos(Q(:,2*i - 1));
        bobs(:,:,i+1) = [xi yi];

    end

    for i = 1:steps

        delete(P);
        clearpoints(h);
        start = max(1,i-L);

        for j = 1:N+1

            addpoints(h, bobs(i,1,j), bobs(i,2,j)); %adds bob points
            %P = [P plot(bobs([start:i],1,j),bobs([start:i],2,j),'Color',[0 j/(N+1) 0])];   %adds bob trajectories for all bobs - LAGGY 

        end

        P = [P plot(bobs([start:i],1,N+1),bobs([start:i],2,N+1),'Color',[0 1 0])];   %adds bob trajectories for end bob only

        drawnow;

    end

end