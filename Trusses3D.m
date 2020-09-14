
clear all
clear screen
close all
clc
%
% This MATLAB code calculates a 2D truss structure
% Code by: Dr. Federica Confalonieri


%Note 1: % is a comment
%Note 2: ; avoids output of the operation
%Note 3: ' after a matrix means transpose
%Note 4: * product between matrices and/or vectors and/or scalars
%Note 5: == means "equal to"
%Note 6: ~= means "different than"

disp('Solution of 3D trusses with the finite element method.')

%Real variable format (alternatives: "format long", "format long e")
format short e

%%Input
%Nodal coordinates
Coor = load('Nodes_coordinates.txt');
Nnodes = size(Coor,1);  %number of nodes
%Connectivity
Connec = load('Elements.txt');
Nelem = size(Connec,1);  %number of elements
%Cross-sectional properties
Sections = load('Sections.txt');
%Boundary conditions
%Nodal forces
Loads = load('Loads.txt');
Nloads = size(Loads,1); %number of loaded nodes
%Restraints (Dirichlet boundary conditions)
bc = load('bc.txt');
Nbc = size(bc,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here ends the input phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Total number of dofs: 3*Nnodes

%Initialization of the global stiffness matrix K 
K=zeros([3*Nnodes,3*Nnodes]);

%Initialization of the right hand side (r.h.s.) vector F
F=zeros([3*Nnodes,1]);

%Initialization of reaction force vector
R=zeros([3*Nnodes,1]);

%Initialization of the axial strains in the trusses
q=zeros([Nelem,1]);

%Initialization of the axial forces in the trusses
Q=zeros([Nelem,1]);
%% ASSEMBLY

%GLOBAL STIFFNESS MATRIX
%Building of matrix K and r.h.s. F
%Loop over elements to assembly K
for ne=1:Nelem
  %element connectivities
  n1 = Connec(ne,1);
  n2 = Connec(ne,2);
  %nodal coordinates
  x1 = Coor(n1,1);
  y1 = Coor(n1,2);
  z1 = Coor(n1,3);
  x2 = Coor(n2,1);
  y2 = Coor(n2,2);
  z2 = Coor(n2,3);
  %calculate element length
  lung=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
  l=sqrt((x2-x1)^2+(y2-y1)^2);
  %calculate sine and cosine of the angle between the horizontal and the
  %truss' axis
  if l==0
      st=0
      ct=0
  else
    st=(y2-y1)/l;
    ct=(x2-x1)/l;
  end 
  
  stq=st^2;
  ctq=ct^2;
  cs=ct*st;
  
  %Computation of the sin and cos of the second considered angle(phi)
  
  sf=(z2-z1)/lung;
  cf=sqrt((x2-x1)^2+(y2-y1)^2)/lung;
  cfq= cf^2;
  sfq=sf^2;
  
   %Local stiffness
  EA = Sections(ne,1)*Sections(ne,2);
  
  %calculate the element stiffness matrix in the global reference frame
    
  kne=(EA/lung)*[cfq*ctq,      cfq*cs,      sf*cf*ct,   -cfq*ctq,   -cfq*cs,    -sf*cf*ct;
                 cfq*cs,       cfq*stq,     sf*cf*st,   -cfq*cs,    -cfq*stq,   -sf*cf*st;
                 sf*cf*ct,     sf*cf*st,    sfq,        -sf*cf*ct,  -sf*cf*st,  -sfq;
                 -cfq*ctq,     -cfq*cs,     -sf*cf*ct,  cfq*ctq,    cfq*cs,     sf*cf*ct;
                 -cfq*cs,      -cfq*stq,    -sf*cf*st,  cfq*cs,     cfq*stq,    sf*cf*st;
                 -sf*cf*ct,    -sf*cf*st,   -sfq,       sf*cf*ct,   sf*cf*st,   sfq;];
 
  
  % look for the correspondence local dofs-global dofs

  global_dofs=[3*n1-2,3*n1-1,3*n1,3*n2-2,3*n2-1,3*n2];
  
  %Loop over element stiffness matrix entries -Does it seem an ASSAMBLY phase?-
  for nl=1:6
    for nm=1:6  
      nr=global_dofs(nl);
      ns=global_dofs(nm);
      K(nr,ns)=K(nr,ns)+kne(nl,nm);
    end
  end
  
end
 
% LOAD VECTOR
%Loop over loaded nodes to assemble F
for ni=1:Nloads
  %loaded node
  n = Loads(ni,1);
  %applied forces
  Px = Loads(ni,2);
  Py = Loads(ni,3);
  Pz = Loads(ni,4)
  %assembly
  dofx=3*n-2;  %global dofs
  dofy=3*n-1;
  dofz=3*n;
  F(dofx)=F(dofx)+Px;
  F(dofy)=F(dofy)+Py;
  F(dofz)=F(dofz)+Pz;
end 

%% BOUNDARY CONDITIONS IMPOSITION
%Modify K because of the restraints
%calculate the norm of diagonal elements in the stiffness matrix
alfa=0;
for ni=1:3*Nnodes
   alfa=alfa+K(ni,ni)^2;
end
alfa=sqrt(alfa)*10^10; %penalty method
%loop over restrained nodes
for ni=1:Nbc
  n = bc(ni,1); %node number
  bcx = bc(ni,2);
  bcy = bc(ni,3);
  bcz = bc(ni,4);
  dofx=3*n-2;  
  dofy=3*n-1;
  dofz=3*n;
  if(bcx==1)
    K(dofx,dofx)=K(dofx,dofx)+alfa;
  end
  if(bcy==1)
    K(dofy,dofy)=K(dofy,dofy)+alfa;
  end
  if(bcz==1)
    K(dofz,dofz)=K(dofz,dofz)+alfa;
  end
end

%% SOLUTION
%Solution of the linear system
u=K\F;
for i=1: length(u)
    if(abs(u(i))< 10^(-10))
        u(i)=0;
    end 
end


%% POST-PROCESSING

%Calculate the reaction forces
%loop over restrained nodes
for ni=1:Nbc
  n=bc(ni,1)
  bcx = bc(ni,2);
  bcy = bc(ni,3);
  bcz = bc(ni,4);
  dofx=3*n-2;  
  dofy=3*n-1;
  dofz=3*n;
  if(bcx==1)
    K(dofx,dofx)=K(dofx,dofx)-alfa; 
    R(dofx)=K(dofx,:)*u-F(dofx);
  end
  if(bcy==1)
    K(dofy,dofy)=K(dofy,dofy)-alfa; 
    R(dofy)=K(dofy,:)*u-F(dofy);
  end
  if(bcz==1)
    K(dofz,dofz)=K(dofz,dofz)-alfa; 
    R(dofz)=K(dofz,:)*u-F(dofz);
  end
end


%Building of the solution for the post-processing
for ne=1:Nelem
  %read element connectivities
  n1=Connec(ne,1);
  n2=Connec(ne,2);
  %read element nodal coordinates
  x1=Coor(n1,1);
  y1=Coor(n1,2);
  z1=Coor(n1,3);
  x2=Coor(n2,1);
  y2=Coor(n2,2);
  z2=Coor(n2,3);
  %calculate element length
  lung=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
  l=sqrt((x2-x1)^2+(y2-y1)^2);
  %calculate sine and cosine of the angle between the horizontal and the
  %truss' axis
  if l==0
      st=0
      ct=0
  else
    st=(y2-y1)/l;
    ct=(x2-x1)/l;
  end 
  
  stq=st^2;
  ctq=ct^2;
  cs=ct*st;
  
  %Computation of the sin and cos of the second considered angle(phi)
  
  sf=(z2-z1)/lung;
  cf=sqrt((x2-x1)^2+(y2-y1)^2)/lung;
  cfq= cf^2;
  sfq=sf^2;
  % look for the correspondence local dofs-global dofs
  global_dofs=[3*n1-2,3*n1-1,3*n1,3*n2-2,3*n2-1,3*n2];
  %build element displacement vector
  uelem=[u(global_dofs(1)),u(global_dofs(2)),u(global_dofs(3)),u(global_dofs(4)),u(global_dofs(5)),u(global_dofs(6))]';
 
  %build rotation matrix, for the 3D case
  
  rot=[cf*ct,   cf*st,  sf,     0,      0,      0;
        0,        0,    0,     cf*ct,  cf*st,  sf];
  
    %build the element displacement vector in the local reference frame
  uloc=rot*uelem;
  %calculate the element axial strain
  qelem=(uloc(2)-uloc(1))/lung;
  
  %calculate element axial force
  EA = Sections(ne,1)*Sections(ne,2);
  Qelem=(EA)*qelem;
  
  %assembly the local axial strains and forces in a global vector 
  q(ne)=qelem;
  Q(ne)=Qelem;
end

%Print results
disp('Nodal displacements u=')
disp(u)
disp('Reaction forces R=')
disp(R)
disp('Axial strains q=')
disp(q)
disp('Axial forces Q=')
disp(Q)


figure(1)
hold on
ampl= 1/max(abs(u));
% Show the mesh 
for ne=1:Nelem
  %read element connectivities
  n1=Connec(ne,1);
  n2=Connec(ne,2);
  %read element nodal coordinates
  x1=Coor(n1,1);
  y1=Coor(n1,2);
  z1=Coor(n1,3);
  x2=Coor(n2,1);
  y2=Coor(n2,2);
  z2=Coor(n2,3);
  % draw the trusses (undeformed shape)
  pA=[x1 x2];
  pB=[y1 y2];
  pC=[z1 z2];
  line(pA,pB,pC,'Color',[0 0 1],'LineWidth',2)
  % draw the trusses (deformed shape)
    % displacement amplification factor
  x1=Coor(n1,1)+ampl*u(3*n1-2);
  y1=Coor(n1,2)+ampl*u(3*n1-1);
  z1=Coor(n1,3)+ampl*u(3*n1);
  x2=Coor(n2,1)+ampl*u(3*n2-2);
  y2=Coor(n2,2)+ampl*u(3*n2-1);
  z2=Coor(n2,3)+ampl*u(3*n2);
  pA=[x1 x2];
  pB=[y1 y2];
  pC=[z1 z2];
  line(pA,pB,pC,'Color',[1 0 0],'LineWidth',2)
end
title('Deformed shape')
hold off

figure(2)
hold on
for ne=1:Nelem
  %read element connectivities
  n1=Connec(ne,1);
  n2=Connec(ne,2);
  %read element nodal coordinates
  x1=Coor(n1,1);
  y1=Coor(n1,2);
  z1=Coor(n1,3);
  x2=Coor(n2,1);
  y2=Coor(n2,2);
  z2=Coor(n2,3);
  % draw the trusses (undeformed shape)
  pA=[x1 x2];
  pB=[y1 y2];
  pC=[z1 z2];
  if (Q(ne)>=0)
    line(pA,pB,pC,'Color',[0 0 1],'LineWidth',2)
  else
    line(pA,pB,pC,'Color',[1 0 0],'LineWidth',2)
  end
  
  
  text((x1+x2)/2,(y1+y2)/2,(z1+z2)/2,num2str(Q(ne)))
end
title('Axial forces')
hold off