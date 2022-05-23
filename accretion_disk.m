npoints = 20;
G = 10;
x_init = 30;

%time parameters
t_final = 10000;
n_steps = 100000;
dt = t_final/n_steps;

%parameters to determine shape of orbit, velocities
M1 = 20;
M2 = 10;
a = 15;
b = 20;
va = sqrt((b/a)*(2*G*(M1+M2))/(a+b));
vb = sqrt((a/b)*(2*G*(M1+M2))/(a+b));

%set field
n_field = 100;
X_field = linspace(-x_init, x_init, n_field);
Y_field = linspace(-x_init, x_init, n_field);
Z_field = linspace(-x_init, x_init, n_field);


%create first body, with initalized position and velocity
M1_radius = 1;
M1_pos = (-M2/(M1+M2))*[a, 0, 0];
M1_vel = (-M2/(M1+M2))*[0, va, 0];
[X1, Y1, Z1] = create_body(npoints, M1_radius);

%create second body, with intialized position and velocity
M2_radius = 10;
M2_density = M2/((4/3)*(M2_radius^3));
M2_pos = (M1/(M1+M2))*[a, 0, 0];
M2_vel = (M1/(M1+M2))*[0, va, 0];
[X2, Y2, Z2] = create_body(npoints, M2_radius);

%Initial value of rate of change of seperation
U = [0,va,0];

particles_x = [];
particles_y = [];
particles_z = [];

particles_u = [];
particles_v = [];
particles_w = [];

for count = 1:t_final
    clf
    X = M2_pos-M1_pos;
    r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);
    
    U = U - dt*G*(M1+M2)*X/(r^3);
    
    M1_vel = -M2/(M1+M2)*U;
    M2_vel = M1/(M1+M2)*U;
    
    M1_pos = M1_pos + dt*M1_vel;
    M2_pos = M2_pos + dt*M2_vel;
    
    %graph spheres centered at current positions
    M1_X = M1_pos(1) + X1;
    M1_Y = M1_pos(2) + Y1;
    M1_Z = M1_pos(3) + Z1;
    
    M2_X = M2_pos(1) + X2;
    M2_Y = M2_pos(2) + Y2;
    M2_Z = M2_pos(3) + Z2;
    
    if (mod(count, 2) == 0)
    surfl(M1_X, M1_Y, M1_Z); hold on
    surfl(M2_X, M2_Y, M2_Z); hold on
    
    phi = zeros(n_field, n_field);
    
    omega = sqrt(G*(M1+M2)/((norm(M2_pos-M1_pos))^3));
    
    %display potential at value in the field
    for r = 1:n_field
        for c = 1:n_field
            phi(r,c) = potential(G, X_field(r), Y_field(c), M1, M2, M1_pos, M2_pos, omega);
        end
    end
    
    L1 = lagrange(G, n_field, M1, M2, M1_pos, M2_pos, omega);
    
    dm = get_mdot(G, M1, M2, M1_pos, M2_pos, M2_radius, L1, omega);
    n_particles = ceil(dm*300);
    
    %intialize pos and vel for discharged particles
    v_particle = -5*(sqrt(M1/M2)/sqrt(norm(M2_pos-M1_pos)/M2_radius))*((M2_pos-M1_pos)/norm(M2_pos-M1_pos));
    
    for n = 1:n_particles
        particles_x(end+1) = M2_pos(1)+randn*M2_radius/10;
        particles_y(end+1) = M2_pos(2)+randn*M2_radius/10;
        particles_z(end+1) = randn*M2_radius/10;
        particles_u(end+1) = v_particle(1);
        particles_v(end+1) = v_particle(2);
        particles_w(end+1) = v_particle(3);
    end
    
    %calculate pos and vel for particles in accretion disk
    for n = 1:length(particles_x)
        
        dx = M1_pos(1)-particles_x(n);
        dy = M1_pos(2)-particles_y(n);
        dz = M1_pos(3)-particles_z(n);
        rr = sqrt(dx.^2 + dy.^2 + dz.^2);
        dxorr3 = dx ./ rr.^3;
        dyorr3 = dy ./ rr.^3;
        dzorr3 = dz ./ rr.^3;
        particles_u(n) = particles_u(n) + dt*G*M1*dxorr3;
        particles_v(n) = particles_v(n) + dt*G*M1*dyorr3;
        particles_w(n) = particles_w(n) + dt*G*M1*dzorr3;
        
        dx = M2_pos(1)-particles_x(n);
        dy = M2_pos(2)-particles_y(n);
        dz = M2_pos(3)-particles_z(n);
        rr = sqrt(dx.^2 + dy.^2 + dz.^2);
        dxorr3 = dx ./ rr.^3;
        dyorr3 = dy ./ rr.^3;
        dzorr3 = dz ./ rr.^3;
        particles_u(n) = particles_u(n) + dt*G*M2*dxorr3;
        particles_v(n) = particles_v(n) + dt*G*M2*dyorr3;
        particles_w(n) = particles_w(n) + dt*G*M2*dzorr3;
        
        particles_x(n) = particles_x(n) + dt*particles_u(n);
        particles_y(n) = particles_y(n) + dt*particles_v(n);
        particles_z(n) = particles_z(n) + dt*particles_w(n);
    end
    
    M2_radius = (.75*(M2-dm)/M2_density)^(1/3); %compute new radius
    M2 = M2-dm;
    M1 = M1+dm;
    [X2, Y2, Z2] = create_body(npoints, M2_radius);
    
    surfl(X_field, Y_field, phi); hold on
    scatter3(particles_x, particles_y, particles_z); hold on
    
    
    xlim([-30, 30])
    ylim([-30, 30])
    zlim([-50, 10])
    
    pause(.1);
    drawnow
    end
end

%create a sphere to represent a body
function [Xcoords, Ycoords, Zcoords] = create_body(npoints, radius)

    Xcoords=zeros(npoints, npoints);
    Ycoords=zeros(npoints, npoints);
    Zcoords=zeros(npoints, npoints);
    
    for r = 1:npoints
        for c = 1:npoints
            Xcoords(r,c) = radius*cos((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Ycoords(r,c) = radius*sin((2*pi*c)/(npoints-1))*sin(pi*r/(npoints-1));
            Zcoords(r,c) = -1*radius*cos((pi*r)/(npoints-1));
        end
    end
end

%calculate gravitational and centrifugal potential
function phi = potential(G, X, Y, M1, M2, M1_pos, M2_pos, omega)

    CM = (M1*M1_pos + M2*M2_pos)/(M1+M2);
    
    phi = -G*M1/(sqrt((X-M1_pos(2))^2+(Y-M1_pos(1))^2)) - G*M2/(sqrt((X-M2_pos(2))^2+(Y-M2_pos(1))^2));
    phi = phi - .5*(norm(cross([0, 0, omega], [X, Y, 0]-CM)))^2;
end

%find lagrange point (saddle point of the potential)
function L1 = lagrange(G, n_field, M1, M2, M1_pos, M2_pos, omega)
    
    r = M2_pos - M1_pos;
    dr = 1/n_field;
    
    phis = ones(n_field-2, 1);
    
     for i = 1:(n_field-2)
         phis(i) = potential(G, M1_pos(1)+i*dr*r(1), M1_pos(2)+i*dr*r(2), M1, M2, M1_pos, M2_pos, omega);
     end
        
    saddle = max(phis);
    ival = find(phis==saddle);
    L1 = M1_pos + ival*dr*r;
end

%calculate rate of mass loss
function dm = get_mdot(G, M1, M2, M1_pos, M2_pos, M2_radius, L1, omega)
    const = .5;
    dm = 0;
    if (norm(L1-M1_pos) < norm(M2_pos-M1_pos)-M2_radius)
        return;
    end
    
    delta_r = norm(M2_pos-L1);
    
    dm = const*M2*(delta_r/M2_radius)^3;
end
