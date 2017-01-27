function [vertexDess] = volume_invariant_3d(Mesh, Para)
% Input arguments:
%   Para.gridRes: length of each grid cell
%   Para.rMin: minimum of the convolution radius
%   Para.rMax: maximum of the convolution radius
%   Para.dessDim: dimension of the descriptor
%   All these values are related to the diagonal length of object bounding
%   box
% This 
SURFELS = mesh_2_surfel(Mesh);
GRID = volumetric_representation(SURFELS, Para);
end

function [SURFELS] = mesh_2_surfel(Mesh)
% Convert a triangular mesh into a surfel representation
p1 = Mesh.vertexPoss(:, Mesh.faceVIds(1,:));
p2 = Mesh.vertexPoss(:, Mesh.faceVIds(2,:));
p3 = Mesh.vertexPoss(:, Mesh.faceVIds(3,:));
% Compute face normals
e12 = p1 - p2;
e13 = p1 - p3;
faceNors = cross(e12, e13);
norms = sqrt(sum(faceNors.*faceNors));
faceNors = faceNors./kron(ones(3,1), norms);
%
numV = size(Mesh.vertexPoss, 2);
numF = size(Mesh.faceVIds, 2);
FVMap = sparse(ones(3,1)*(1:numF), Mesh.faceVIds, ones(3, numF), numF, numV);
%
vertexNors = double(faceNors)*FVMap;
vertexNors = single(vertexNors);
SURFELS = [Mesh.vertexPoss; vertexNors];
end

function [GRID] = volumetric_representation(SURFELS, Para)
% Compute the volumetric representation of a mesh.
% We first convert the mesh into a point cloud, then for every center of
% the grid cell, we compute its distance to the point cloud. We will use
% that information to build the grid representation.
vertexPoss = SURFELS(1:3, :);
vertexNors = SURFELS(4:6, :);
lowerCorner = min(vertexPoss')';
upperCorner = max(vertexPoss')';
boxDim = double(upperCorner - lowerCorner);
diagonalLength = norm(boxDim);
gridRes = double(Para.gridRes*diagonalLength);
rMin = double(Para.rMin*diagonalLength);
rMax = double(Para.rMax*diagonalLength);
dessDim = Para.dessDim;
dimX = floor(boxDim(1)/gridRes) + 2;
dimY = floor(boxDim(2)/gridRes) + 2;
dimZ = floor(boxDim(3)/gridRes) + 2;
center = double(lowerCorner + upperCorner)/2;
GRID.lowerCorner = [center(1) - dimX*gridRes/2;
    center(2) - dimY*gridRes/2;
    center(3) - dimZ*gridRes/2];
GRID.dimX = dimX;
GRID.dimY = dimY;
GRID.dimZ = dimZ;
% Locations of cell centers
cellCenters = zeros(3, dimX*dimY*dimZ);
cellCenters(1,:) = (GRID.lowerCorner(1)+0.5*gridRes)*ones(1, dimX*dimY*dimZ)...
    + kron(ones(1,dimY*dimZ), 0:(dimX-1))*gridRes;
cellCenters(2,:) = (GRID.lowerCorner(2)+0.5*gridRes)*ones(1, dimX*dimY*dimZ)...
    + kron(kron(ones(1,dimZ), 0:(dimY-1)), ones(1,dimX))*gridRes;
cellCenters(3,:) = (GRID.lowerCorner(3)+0.5*gridRes)*ones(1, dimX*dimY*dimZ)...
    + kron(0:(dimZ-1), ones(1,dimX*dimY))*gridRes;
cellCenters = single(cellCenters);
IDX = knnsearch(vertexPoss', cellCenters')';
pc_dis = sum((cellCenters - vertexPoss(:, IDX)).*vertexNors(:, IDX));
pc_dis = pc_dis/gridRes;
ids_0 = find(pc_dis > 1);
ids_1 = find(pc_dis <= 1 & pc_dis >= -1);
ids_2 = find(pc_ids < -1);
pc_dis(ids_0) = 0;
pc_dis(ids_2) = 1;
pc_dis(ids_1) = (1 - pc_dis(ids_1))/2;
GRID.data = reshape(pc_dis, [dimX, dimY, dimZ]);
GRID.rMin = rMin/gridRes;
GRID.rMax = rMax/gridRes;
GRID.dessDim = dessDim;
end

function [grid3d_k] = convolution_kernel(GRID)
% Compute the convolution kernel of a sphere with radius on a grid mesh
% with resolution gridRes

end

function [vertexDess] = interpolate_descriptors(Mesh, GRID)
% Interpolate the volumetric descriptor to obtain descriptors at mesh
% vertices
end