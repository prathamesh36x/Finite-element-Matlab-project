N = 5;
conn = [[1:N-1;2:N] [N+1:2*N-2;N+2:2*N-1] [1:N-1;N+1:2*N-1] [N+1:2*N-1;2:N]];
conn = [conn (2*N-1)+[[1:N-1;2:N] [N+1:2*N-2;N+2:2*N-1] [1:N-1;N+1:2*N-1] [N+1:2*N-1;2:N]]];
conn = [conn [[1:N N+1:2*N-1];(2*N-1)+[1:N N+1:2*N-1]]];
conn = [conn [1:N-1;3*N-1+(1:N-1)] [2:N;2*N:3*N-2]]; %[2*N:3*N-2;N+1:2*N-1] 
x = [0:N-1 (0:N-2)+0.5; 0*(0:N-1) 1+(0:N-2)*0]; %meters
x = [[x;x(1,:)*0] [x;x(1,:)*0+2]]; x = [x(1,:);x(3,:);-x(2,:)]; %swith y and z dir

subplot(2,1,1)
for c = conn
    xe = x(:,c);
    plot3(xe(1,:),xe(2,:),xe(3,:),'k-')
    hold on
end


axis equal

E = 5.7*10^9; %5.7GPa Youngs Modulus carbon fiber
dia = 0.2; %mm == 10mm radius
A = pi*(dia/2)^2; %cross-sectional area
ro = 1.6e3; %kg/m^3
K = zeros(length(x)*3); %stiffness Matrix
f = zeros(length(x)*3,1); %forces
Ke2 = E * A * [1, -1; -1, 1];

for c = conn
    xe = x(:,c);
    dx = xe(:,2) - xe(:,1);
    Re = [dx', 0, 0, 0;  0, 0, 0, dx'] / norm(dx);     
    sctr(1:3:6) = 3*c-2;
    sctr(2:3:6) = 3*c-1;
    sctr(3:3:6) = 3*c;
    K(sctr,sctr) = K(sctr,sctr) + Re' * Ke2 / norm(dx) * Re;
end

all_nodes = 1:length(x);
vs = 6; %m/s crosswind speed 

% Adding force of 1000 Newton in the negative z direction to top nodes
force_magnitude = 1000000;
top_nodes = find(x(3,:) == max(x(3,:)));
f(3*top_nodes) = -force_magnitude;

% Fixing Degrees of freedom
fixed = [find(x(1,:)==0) find(x(1,:)==0.5) find(x(1,:)==N-1.5) find(x(1,:)==(N-1))];
fixed_dof = [3*fixed-2 3*fixed-1 3*fixed];
K(fixed_dof,:) = 0;
K(fixed_dof,fixed_dof) = eye(length(fixed_dof));
f(fixed_dof) = 0;
plot3(x(1,fixed),x(2,fixed),x(3,fixed),'ko') % Indicate fixed Locations

quiver3(x(1,all_nodes),x(2,all_nodes),x(3,all_nodes), ...
    f(3*all_nodes-2)',0*f(3*all_nodes-1)',f(3*all_nodes)',0.5,'color','red')
quiver3(x(1,all_nodes),x(2,all_nodes),x(3,all_nodes), ...
    0*f(3*all_nodes-2)',f(3*all_nodes-1)',0*f(3*all_nodes)',0.5,'color','blue')

[X,Y] = meshgrid(-1:1:N,-1:1:3);
Z = -sign(X) + sign(X-N+1); Z(Z==-1) = 0;
surf(X,Y,Z,abs(gradient(Z))); alpha 0.5; camproj('perspective'); view([0 20])
title('Truss Structure and Force Visualization - Red = Gravity & Blue = Applied Force'); hold off

% Solving for displacements
d = K\f;

% Update node positions with displacements
xn = x + [d(1:3:end), d(2:3:end), d(3:3:end)]';

% Plotting deformed structure
subplot(2,1,2); 
AF = 0.025/max(abs(d)); % Deformation amplification factor
xn = x + AF*[d(1:3:end), d(2:3:end), d(3:3:end)]';
for i = 1:length(conn)
    c = conn(:,i);
    xe = xn(:,c);
    plot3(xe(1,:),xe(2,:),xe(3,:),'.-','markersize',15,'MarkerEdgeColor','black',...
        'linewidth',2,'Color','blue');  hold on
end
hold off

grid minor
surface(X,Y,Z); alpha 0.5; camproj('perspective'); view([-13 11])
axis([0 N-1 -.1 2.1 -2.1 .1]); title('Deformation (Exaggerated)');
axis equal
disp(d)






clear all;
close all;
clc;

% Material properties
E = 210e9; % Young's modulus

% Truss geometry
nodes = [
    0, 0, 0;
    5, 0, 0;
    10, 0, 0;
    0, 5, 0;
    5, 5, 0;
    10, 5, 0;
    0, 0, 5;
    5, 0, 5;
    10, 0, 5;
];

% Truss connectivity (8-noded truss elements)
elements = [
    1, 2;
    2, 3;
    4, 5;
    5, 6;
    7, 8;
    8, 9;
    1, 4;
    2, 5;
    3, 6;
    4, 7;
    5, 8;
    6, 9;
];

% Applied loads
applied_forces = zeros(size(nodes));
applied_forces(6, 3) = -1000; % Applied force at node 6 in the z-direction

% Gauss quadrature points and weights for 2-point quadrature
gauss_points = [-sqrt(1/3), sqrt(1/3)];
gauss_weights = [1, 1];

% Initialize global stiffness matrix and load vector
num_nodes = size(nodes, 1);
num_dofs_per_node = 3;
K_global = zeros(num_nodes * num_dofs_per_node);
F_global = zeros(num_nodes * num_dofs_per_node, 1);

% Loop over elements to assemble global stiffness matrix and load vector
for i = 1:size(elements, 1)
    element_nodes = elements(i, :);
    node_coords = nodes(element_nodes, :);
    
    % Compute element stiffness matrix using Gauss quadrature
    k_element = zeros(num_dofs_per_node * 2);
    for j = 1:length(gauss_points)
        xi = gauss_points(j);
        [N, dN_dxi] = shape_functions(xi);
        
        % Jacobian matrix
        J = node_coords' * dN_dxi';
        
        % Derivative of shape functions with respect to physical coordinates
        dN_dx = dN_dxi / J;
        
        % Element stiffness matrix for one Gauss point
        B = compute_strain_displacement_matrix(dN_dx);
        k_element = k_element + gauss_weights(j) * B' * E * A * B * det(J);
    end
    
    % Assemble element into global stiffness matrix and load vector
    indices_local = [num_dofs_per_node * element_nodes(1) - (num_dofs_per_node - 1):num_dofs_per_node * element_nodes(1), ...
                     num_dofs_per_node * element_nodes(2) - (num_dofs_per_node - 1):num_dofs_per_node * element_nodes(2)];

    indices_global = indices_local + (i - 1) * num_dofs_per_node;
    
    K_global(indices_global, indices_global) = K_global(indices_global, indices_global) + k_element;
    
    % Compute element load vector (assuming constant distributed load)
    f_element = zeros(num_dofs_per_node * 2, 1);
    f_element(end - 2:end) = [1; 1; 1] * applied_forces(element_nodes(2), 3) / length(gauss_points);
    
    F_global(indices_global) = F_global(indices_global) + f_element;
end

% Apply boundary conditions (fix nodes at the base)
fixed_nodes = [1, 2, 3, 7, 8, 9]; % Fix nodes at the base in all directions
fixed_dofs = reshape(repmat(fixed_nodes, 1, num_dofs_per_node), 1, []);
K_global(fixed_dofs, :) = [];
K_global(:, fixed_dofs) = [];
F_global(fixed_dofs) = [];

% Solve for displacements
displacements = K_global \ F_global;

% Assemble complete displacement vector
full_displacements = zeros(num_nodes * num_dofs_per_node, 1);
full_displacements(setdiff(1:length(full_displacements), fixed_dofs)) = displacements;

% Display results
disp('Nodal Displacements:');
disp(reshape(full_displacements, num_dofs_per_node, []).');

