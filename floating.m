clear all
clc
stlFile = 'cad/Part1.STL';

%% params
E = 8e10*1e-6;     % Young's modulus in MPa
nu = 0.3;     % Poisson's ratio (nondimensional)
rho = 2700*1e-12;    % Mass density in t/mm^3
gravity = [0 0 -9.80665];

origins = [0 0 1567/2];

numFrames = size(origins,1);

femodel = createpde('structural','modal-solid');
importGeometry(femodel,stlFile);

structuralProperties(femodel, ...
 'YoungsModulus',E, ...
 'PoissonsRatio',nu, ...
 'MassDensity',rho);

generateMesh(femodel, ...
 'GeometricOrder','quadratic','Hmax',1e2,'Hmin',1e2);

faceIDs = [3];    % List in the same order as the interface frame origins

structuralBC(femodel, ...
 'Face',faceIDs(1), ...
 'Constraint','multipoint', ...
 'Reference',origins(1,:));

rom = reduce(femodel,'FrequencyRange',[0 2]*(2*pi));
TransformationMatrix=...
rom.MPCTransformationMatrix*rom.TransformationMatrix; % refer to  which pde.ReducedStructuralModel

%% compare with ansys
% Test modal frequency
% fix - free 1.6708 Hz
assert(abs(sqrt(rom.K(end,end))/2/pi-1.6708)<0.2)


volumes = zeros(1,size(femodel.Mesh.Elements,2));
for i= 1:size(femodel.Mesh.Elements,2)
    volumes(i) = volume(femodel.Mesh,i)/volume(femodel.Mesh);
end

weights= zeros(1,size(femodel.Mesh.Nodes,2));
for i= 1:size(femodel.Mesh.Nodes,2)
    weights(i)=sum(volumes(sum(femodel.Mesh.Elements==i)==1))/10;
end
m_i = volume(femodel.Mesh)*femodel.MaterialProperties.MaterialAssignments.MassDensity*1e3;

digitsold = digits(7);
% Generalized coordinates:
% q1, q2, q3: 3D position (x, y, z)
% q4, q5, q6, q7: Quaternion (q, qx, qy, qz)
% q8: Modal coordinate (flexible deformation)
syms q1 q2 q3 q4 q5 q6 q7 q8 
syms dq1 dq2 dq3 w1 w2 w3 dq8 

q = [q1, q2, q3,... 
    q4, q5, q6, q7,... 
    q8].';
dq = [dq1, dq2, dq3,... 
    w1, w2, w3,... 
    dq8].';

x=[0, 0, 0, 0, 0, 0, q(end)].'; 

T7 = 0;
V7 = 1/2*x.'*1e3*sym(rom.K,'d')*x;

AttachmentDoF = rom.RetainedDoF;
CondensedDOF = true(3*size(femodel.Mesh.Nodes,2), 1);
CondensedDOF(AttachmentDoF) = false;
CondensedDoF = find(CondensedDOF);
DOF = [AttachmentDoF, CondensedDoF'];

disp=sym(TransformationMatrix,'d')*x;

rot=setOrientation([q4, q5, q6, q7]);
d= [q1, q2, q3];
r= d + [0, 0, -1567/2];
omega = [w1, w2, w3];
v= [dq1, dq2, dq3];

for i= 1:size(femodel.Mesh.Nodes,2)
    s_i = sym(femodel.Mesh.Nodes(:,i).','d')-sym(origins(1,:),'d');
    u_i = [disp(DOF==i),disp(DOF==i+size(femodel.Mesh.Nodes,2)),disp(DOF==i+2*size(femodel.Mesh.Nodes,2))];
    r_i = d + (s_i*1e-3+u_i);
    % Elastic velocity: partial derivative of u_i w.r.t q8 * dq8
    v_elastic_local = jacobian(u_i, q(8)) * dq(end);
    v_i = v + ...                          % Translational Velocity
          cross(omega, r) + ...            % Rotational Velocity (w x r)
          v_elastic_local.'*rot.';       % Elastic Velocity (Rotated to global)
    
    T7 = T7 + 1/2*sym(m_i*weights(i),'d')*(v_i*v_i.');
    mg_i = sym(m_i*weights(i),'d')*gravity*rot.';
    V7 = V7 - mg_i*(d + s_i*1e-3).';
end

% Compute EoM
T = T7;
V = V7;

% use the chain rule: d/dt(f) = (df/dq)*dq_actual
% Kinematic Equation: q_dot = 0.5 * E(q)' * w
% Using the E matrix definition 
    E = [-q(2),  q(1),  q(4), -q(3);
         -q(3), -q(4),  q(1),  q(2);
         -q(4),  q(3), -q(2),  q(1)];
% IMPORTANT: dq_actual must use the quaternion kinematic relation for the rotational part
dq_actual = [dq(1:3); % v
             0.5 * E' * dq(4:6); % q_dot (4x1)
             dq(7)]; % dq8
% Projecting the 8-dim gradient back to the 7-dim quasi-velocity space:
% For translation and modal: identity mapping
% For rotation: Use the 0.5 * E mapping
W_mat = blkdiag(eye(3), 0.5 * E, 1); % This is a 8x7 mapping matrix

M = jacobian(jacobian(T, dq), dq);
Cq = ((jacobian(jacobian(T,dq),q))*dq_actual-W_mat*jacobian(T,q).');
K = (W_mat*jacobian(V,q).');

Z = [q ; dq];

M_IMM = [eye(length(q)) zeros(length(q),size(M,1));
       zeros(size(M,1),length(q)) M];

F_in = [zeros(8,1);1;zeros(6,1)];

F_sys= -Cq-K;

gbar=[zeros(length(dq),1);F_in];
fbar=[dq;F_sys];

%% Define control Affine system \dot{x}=f(x)+g(x)u
f=pinv(M_IMM)*fbar;

g=pinv(M_IMM)*gbar;

% ffunc = matlabFunction(f,'File','ffunc',"Vars",{Z});
% gfunc = matlabFunction(g,'File','gfunc',"Vars",{Z});