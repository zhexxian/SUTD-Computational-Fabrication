% MATLAB file for finite element solution of 2D plane stress problem,
% thin plate under traction force
% 
% Singapore University of Technology and Design
% Course: ISTD 01.110, Computational Fabrication, 2017
%
% Name: Zhang Zhexian
% Student ID: 1001214
%

function fem_2d_plate() 

%--- Input parameters

% Material parameters 
E  = 10e7;          % Pa = N/m²  (Young’s modulus)
nu = 0.3;          % (Poisson’s ratio)
thickness = 0.1;    % m

% Traction force
P = 1e7;    % N/m (on the boundary)

% Mesh dimensions and size
Lx = 10;     % m (length)
Ly = 5;     % m (width)
numberElementsX = 40;
numberElementsY = 20;

% Figure for plots
figure; hold on;
axis equal;
axis([-2 12 -2 7]);


%--- 1. Mesh generation
%       Create rectangular mesh

% * Nodes
nodeSpacingX = Lx/numberElementsX;    %(ZX)  % Spacing of nodes in x-direction
nodeSpacingY = Ly/numberElementsY;    %(ZX)  % Spacing of nodes in y-direction
numberNodes = (numberElementsX + 1) * (numberElementsY + 1);
nodeCoordinates = zeros(numberNodes,2);
    % nN x 2 matrix containing the (x,y) coordinates of the nN nodes
nn = 1;
for ny = 0:numberElementsY
    for nx = 0:numberElementsX
        nodeCoordinates(nn, :) = [nx*nodeSpacingX, ny*nodeSpacingY];
        nn = nn + 1;
    end;
end;
GDof = 2*numberNodes;   
    % No. of global degrees of freedom, 
    % 2 degrees of freedom for every node (u and v)

% * Elements
numberElements = numberElementsX * numberElementsY;
elementNodes = zeros(numberElements, 4);
    % el x 4 matrix containing the 4 node indices of each of the nE elements
nn = 1;
for ey = 1:numberElementsY
    for ex = 1:numberElementsX
        elementNodes(nn, :) = ...
            [ (ey-1)*(numberElementsX + 1) + ex, ...
              (ey-1)*(numberElementsX + 1) + (ex+1), ...
              ey*(numberElementsX + 1) + (ex+1), ...
              ey*(numberElementsX + 1) + ex];
        nn = nn + 1;
    end;
end;

plotMesh(nodeCoordinates, elementNodes);


%--- 2. Discretization of elements
%       see 'shapeFunctionQ4' and 'Jacobian'


%--- 3.+4. Evaluation of weak form and assembly

% Constitutive matrix
C = [(E/(1-nu^2)), ((nu*E)/(1-nu^2)), 0;...
     ((nu*E)/(1-nu^2)),(E/(1-nu^2)), 0; ...
     0, 0, (E/(2*(1+nu)))];  %(ZX)   % Assign the 3x3 constitutive matrix of the 2D plane stress problem

% Assembly of global stiffness matrix

stiffMat = assembleStiffMat2D(elementNodes, nodeCoordinates, C, thickness);

% * Remark:
% The global stiffness matrix (stiffMat), force vector (forceVec) and
% displacement vector (ufull) have size GDof=2*numberNodes and the DOFs are
% ordered such that the [1:numberNodes] values relate to the
% u-displacements of the nodes and the [numberNodes+1:GDof] values relate
% to the v-displacements of the nodes.

% Global force vector (traction force P applied at X==Lx in X-direction)
% * This is a very simple way of assigning the boundary traction forces
%   without looping over all elements, which works if the traction force is
%   constant over the boundary edge.
forceVec = zeros(GDof,1);   
rightBoundary = find(nodeCoordinates(:,1)==Lx);
forceVec(rightBoundary) = P * thickness * nodeSpacingY;
forceVec(rightBoundary(1)) = P * thickness * nodeSpacingY/2;
forceVec(rightBoundary(end)) = P * thickness * nodeSpacingY/2;



% Essential boundary conditions
fixedNodeX = find(nodeCoordinates(:,1)==0); % fixed U where X==0
fixedNodeY = find(nodeCoordinates(:,2)==0); % fixed V where Y==0
prescribedDofs = [fixedNodeX; fixedNodeY+numberNodes];
activeDofs = setdiff(1:GDof, prescribedDofs); 
    % Remove prescribedDofs from the array of all DOFs [1:GDof]


%--- 5. Solution of linear system

% Apply B.C. to global stiffness matrix  
K = stiffMat(activeDofs,activeDofs);
b = forceVec(activeDofs);

% Solve
u = K\b;

% Full displacement vector for all DOFs and nodal displacements
ufull = zeros(GDof,1);
ufull(activeDofs) = u;
nodalDispl = reshape(ufull, numberNodes, 2);

%--- 6. Post-processing

stress = stresses2D(elementNodes, nodeCoordinates, nodalDispl, C);

plotMesh(nodeCoordinates, elementNodes, nodalDispl, stress);
plotMesh(nodeCoordinates, elementNodes);
colorbar;
colormap 'jet';
title('Sigma XX stress (on deformed shape)');

% Energy
%   compute strain energy here (can be done in 1 line of code)
strainEnergyX = 0.5*thickness*transpose(nodalDispl(:,1))*nodalDispl(:,1);
strainEnergyY = 0.5*thickness*transpose(nodalDispl(:,2))*nodalDispl(:,2);
strainEnergy = strainEnergyX + strainEnergyY;
output = ['Strain energy: ',num2str(strainEnergy)];
disp(output);
% Poisson's ratio
%   compute Poisson's ratio based on nodal displacement of top right corner
%   (can be done in 1 line of code)
poissonRatio = - (nodalDispl(end,2)/Ly)/(nodalDispl(end,1)/Lx);
output = ['Poissons ratio: ',num2str(poissonRatio)];
disp(output);
end

% -------------------------------------------------------------------------
% Additional functions that require manual modifications 
% -------------------------------------------------------------------------

function [stiffMat] = assembleStiffMat2D( ...
    elementNodes, nodeCoordinates, C, thickness)
% Compute and assemble stiffness matrix for plane stress Q4 elements
    
    % Initialization
    numberElements = size(elementNodes, 1);
    numberNodes = size(nodeCoordinates, 1);
    GDof = 2*numberNodes;
    stiffMat = zeros(GDof, GDof);

    % 2 by 2 quadrature points
    [gaussWeights, gaussLocations] = gaussQuadrature('complete');
    
    % Loop over elements
    for e = 1:numberElements
        
        indices = elementNodes(e,:);    % indices of nodes of the element
        nind = length(indices);         % no. of element nodes (=4)   
        elementDofs = [indices, indices+numberNodes];   
            % DOFs of the element in the global stiffness matrix, 
            % first part relates to u- and second part to v-displacements    
        ndof = length(elementDofs);     % no. of element DOFs (=2*4)
        % Initialize element stiffness matrix
        elemStiffMat = zeros(ndof,ndof);   %(ZX)    
        % Loop over Gauss points
        for q = 1:size(gaussWeights,1)
            
            % Gauss point location in [-1,1]x[-1,1]
            GaussPoint = gaussLocations(q,:);
            xi = GaussPoint(1);
            eta = GaussPoint(2);
            
            % evaluate shape functions and derivatives at Gauss point
            [shapeFunctions, naturalDerivatives] = shapeFunctionQ4(xi, eta);
            
            % Jacobian matrix and conversion of derivatives from local to
            % global
            [JacobianMat, invJacobian, XYderivatives] = ...
                Jacobian(nodeCoordinates(indices,:), naturalDerivatives);
            
            % B matrix (strain-displacement matrix)
            B = zeros(3,ndof);
            transposeXYderivatives = transpose(XYderivatives);
            B(1,1:nind) = transposeXYderivatives(1,1:nind);       %(ZX)   % Assign the components of the B-matrix
            B(2,nind+1:ndof) = transposeXYderivatives(2,1:nind);   %(ZX)
            B(3,1:nind) = transposeXYderivatives(2,1:nind);         %(ZX)
            B(3,nind+1:ndof) = transposeXYderivatives(1,1:nind);   %(ZX)
            
            % Add evaluation at Gauss point to element stiffness matrix
            elemStiffMat = elemStiffMat + ...
                    thickness*transpose(B)*C*B*det(JacobianMat); %(ZX)    % Compute the integrand of the element stiffness matrix 
                    % integral(thickness*transpose(B)*C*B*det(JacobianMat),-1,1);
        end;
        
        % Add element stiffness matrix to global stiffness matrix
        stiffMat(elementDofs, elementDofs) = ...
                stiffMat(elementDofs, elementDofs) + elemStiffMat;
    end;
    
end

% -------------------------------------------------------------------------

function [shapeFunctions, naturalDerivatives] = shapeFunctionQ4(xi, eta)
% Shape functions and derivatives for Q4 elements
% shapeFunctions : Shape functions
% naturalDerivatives : derivatives w.r.t. xi and eta
% xi, eta : natural coordinates (-1 ... +1)

    shapeFunctions = ...
        1/4 * [ (1-xi)*(1-eta);
                (1+xi)*(1-eta);
                (1+xi)*(1+eta);
                (1-xi)*(1+eta)]; %(ZX)   % evaluate the shape functions N_i^e(xi,eta) into a 4x1 vector
            
    naturalDerivatives = ...
        1/4 * [ -(1-eta), -(1-xi);
                  1-eta,  -(1+xi);
                  1+eta,    1+xi;
                -(1+eta),   1-xi ];

end 

% -------------------------------------------------------------------------

function [JacobianMatrix, invJacobian, XYDerivatives] = ...
    Jacobian(nodeCoordinates, naturalDerivatives)
% Jacobian matrix, inverse Jacobian and natural derivatives for an element
% JacobianMatrix : Jacobian matrix
% invJacobian : inverse of Jacobian Matrix
% XYDerivatives : shape function derivatives w.r.t. x and y
% nodeCoordinates : nodal coordinate matrix at element level
% naturalDerivatives : shape function derivatives w.r.t. xi and eta

    JacobianMatrix = transpose(nodeCoordinates)*naturalDerivatives; %(ZX)  % Use nodal coordinates and natural derivatives to compute the Jacobian matrix J
    invJacobian = inv(JacobianMatrix);
    XYDerivatives = naturalDerivatives * invJacobian;

end 

% -------------------------------------------------------------------------

function stress = stresses2D(elementNodes, nodeCoordinates, nodalDispl, C)

    % Initialization
    numberElements = size(elementNodes, 1);
    numberNodes = size(nodeCoordinates, 1);
    GDof = 2*numberNodes;
    displVec = reshape(nodalDispl, GDof, 1);
    
    % stresses at nodes
    stress = zeros(numberElements,size(elementNodes,2),3);
    
    % Evaluate stresses at element corner nodes
    stressPoints = [-1 -1; 1 -1; 1 1; -1 1];
    
    % Loop over elements
    for e = 1:numberElements
        
        indices = elementNodes(e,:);
        nind = length(indices);
        elementDofs = [indices, indices+numberNodes];
        ndof = length(elementDofs);
        
       
        % Loop over evaluation points
        for q = 1:size(stressPoints,1) 
            
            pt = stressPoints(q,:);           
            xi = pt(1);
            eta = pt(2);

            % shape functions and derivatives
            [shapeFunctions, naturalDerivatives] = shapeFunctionQ4(xi,eta);
            
            % Jacobian etc.
            [JacobianMat, invJacobian, XYderivatives] = ...
                Jacobian(nodeCoordinates(indices,:), naturalDerivatives);
            
            % B matrix (strain-displacement matrix)
            B = zeros(3,ndof);
            transposeXYderivatives = transpose(XYderivatives);
            B(1,1:nind) = transposeXYderivatives(1,1:nind);        %(ZX)   % Assign the components of the B-matrix
            B(2,nind+1:ndof) = transposeXYderivatives(2,1:nind);   %(ZX)   % - same as above in 'assembleStiffMat2D'
            B(3,1:nind) = transposeXYderivatives(2,1:nind);         %(ZX)
            B(3,nind+1:ndof) = transposeXYderivatives(1,1:nind);   %(ZX)           
            
            % strain and stress
            strain = B * displVec(elementDofs);
            stress(e,q,:) = C*strain;      %(ZX)   % Compute stress vector from strain vector
        end;
    end;    
    
end

% -------------------------------------------------------------------------
% Below functions do not have any placeholders and  
% do not require any modifications 
% -------------------------------------------------------------------------

function plotMesh(nodes, elements, def, stress)
% Plot 2D mesh
% nodes : matrix of nodal coordinates
% elements : matrix of element node numbers
% def : matrix of nodal deformations

    if (nargin < 3)
        nodesDef = nodes;
        color = [0.4 0.4 0.4];
    else
        nodesDef = nodes + def;
        color = [0 0 0];
    end;

    for e = 1:size(elements,1)
        
        if (nargin > 3)
            XX = nodesDef(elements(e,:),1);
            XX = XX([1 2; 4 3]);
            YY = nodesDef(elements(e,:),2);
            YY = YY([1 2; 4 3]);
            SS = stress(e,:,1);
            SS = SS([1 2; 4 3]);
            pcolor(XX, YY, SS);
        else
            elnodes = [elements(e,:), elements(e,1)];
            plot(nodesDef(elnodes, 1), nodesDef(elnodes, 2), 'Color', color);
        end;
    end;  

end

% -------------------------------------------------------------------------

function [weights, locations] = gaussQuadrature(option)
% Gauss quadrature for Q4 elements
% option ‘complete’ (2x2)
% option ‘reduced’ (1x1)
% locations : Gauss point locations
% weights : Gauss point weights

    switch option
        case 'complete'
            locations = ...
                [ -0.577350269189626, -0.577350269189626;
                   0.577350269189626, -0.577350269189626;
                   0.577350269189626,  0.577350269189626;
                  -0.577350269189626,  0.577350269189626];
            weights = [ 1;1;1;1];
        case 'reduced'
            locations = [0 0];
            weights = [4];
    end;

end



