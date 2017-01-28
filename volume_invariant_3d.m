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
fprintf('Constructing the volumetric representation...\n');
GRID = volumetric_representation(SURFELS, Para);
fprintf('Performing the FFT Transform of the data structure...\n');
shape_freq = fftn(GRID.data);

vertexDess = single(zeros(Para.dessDim, size(SURFELS, 2)));

for id = 1 : Para.dessDim
    t = (id-1)/(GRID.dessDim -1);
    radius = GRID.rMin*(1-t) + GRID.rMax*t;
    grid_k = convolution_kernel(...
        GRID.dimX, GRID.dimY, GRID.dimZ, radius);
    kernel_freq = fftn(grid_k);
    dess_freq = shape_freq.*kernel_freq;
    dess_vol = ifftn(dess_freq);
    vertexDess_per_radius = interpolate_descriptors(Mesh, GRID, dess_vol);
    vertexDess(id, :) = vertexDess_per_radius;
    fprintf('Finished descriptor_%d\n', id);
end
% Return the vertex descriptor
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
% Remove empty faces
ids = find(norms > 1e-16);
faceNors(:, ids) = faceNors(:,ids)./kron(ones(3,1), norms(ids));
%
numV = size(Mesh.vertexPoss, 2);
numF = size(Mesh.faceVIds, 2);
FVMap = sparse(ones(3,1)*(1:numF), Mesh.faceVIds, ones(3, numF), numF, numV);
%
vertexNors = double(faceNors)*FVMap;
norms = sqrt(sum(vertexNors.*vertexNors));
ids = find(norms > 1e-16);
vertexNors(:, ids) = vertexNors(:,ids)./kron(ones(3,1), norms(ids));
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
rMin = double(Para.rMin*diagonalLength)/gridRes;
rMax = double(Para.rMax*diagonalLength)/gridRes;
dessDim = Para.dessDim;
dimX = floor(boxDim(1)/gridRes) + 4 + 2*floor(rMax);
dimY = floor(boxDim(2)/gridRes) + 4 + 2*floor(rMax);
dimZ = floor(boxDim(3)/gridRes) + 4 + 2*floor(rMax);
center = double(lowerCorner + upperCorner)/2;
GRID.lowerCorner = [center(1) - dimX*gridRes/2;
    center(2) - dimY*gridRes/2;
    center(3) - dimZ*gridRes/2];
GRID.dimX = dimX;
GRID.dimY = dimY;
GRID.dimZ = dimZ;
GRID.gridRes = gridRes;
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
ids_2 = find(pc_dis < -1);
pc_dis(ids_0) = 0;
pc_dis(ids_2) = 1;
pc_dis(ids_1) = (1 - pc_dis(ids_1))/2;
GRID.data = reshape(pc_dis, [dimX, dimY, dimZ]);
GRID.rMin = rMin;
GRID.rMax = rMax;
GRID.dessDim = dessDim;
end

function [grid3d_k] = convolution_kernel(dimX, dimY, dimZ, radius)
% Compute the convolution kernel of a sphere with radius on a grid mesh
% with resolution gridRes
cellCenters = zeros(3, dimX*dimY*dimZ);
cellCenters(1,:) = -(dimX/2-0.5)*ones(1, dimX*dimY*dimZ)...
    + kron(ones(1,dimY*dimZ), 0:(dimX-1));
cellCenters(2,:) = -(dimY/2-0.5)*ones(1, dimX*dimY*dimZ)...
    + kron(kron(ones(1,dimZ), 0:(dimY-1)), ones(1,dimX));
cellCenters(3,:) = -(dimZ/2-0.5)*ones(1, dimX*dimY*dimZ)...
    + kron(0:(dimZ-1), ones(1,dimX*dimY));
cellCenters = single(cellCenters);
pc_dis = sqrt(sum(cellCenters.*cellCenters))-radius;
ids_0 = find(pc_dis > 1);
ids_1 = find(pc_dis <= 1 & pc_dis >= -1);
ids_2 = find(pc_dis < -1);
pc_dis(ids_0) = 0;
pc_dis(ids_2) = 1;
pc_dis(ids_1) = (1 - pc_dis(ids_1))/2;
pc_dis = pc_dis/sum(pc_dis);
grid3d_k = reshape(pc_dis, [dimX, dimY, dimZ]);
end

function [vertexDess] = interpolate_descriptors(Mesh, GRID, dess)
% Interpolate the volumetric descriptor to obtain descriptors at mesh
% vertices
coordsX = (Mesh.vertexPoss(1,:)-GRID.lowerCorner(1))/GRID.gridRes - 0.5;
coordsY = (Mesh.vertexPoss(2,:)-GRID.lowerCorner(2))/GRID.gridRes - 0.5;
coordsZ = (Mesh.vertexPoss(3,:)-GRID.lowerCorner(3))/GRID.gridRes - 0.5;
%
idsX = floor(coordsX);
tX = coordsX -  idsX;
idsY = floor(coordsY);
tY = coordsY -  idsY;
idsZ = floor(coordsZ);
tZ = coordsZ -  idsZ;
%
ids000 = (idsX + 1) + (idsY + idsZ*GRID.dimY)*GRID.dimX;
ids100 = (idsX + 2) + (idsY + idsZ*GRID.dimY)*GRID.dimX;
ids010 = (idsX + 1) + ((idsY+1) + idsZ*GRID.dimY)*GRID.dimX;
ids110 = (idsX + 2) + ((idsY+1) + idsZ*GRID.dimY)*GRID.dimX;
%
ids001 = (idsX + 1) + (idsY + (idsZ+1)*GRID.dimY)*GRID.dimX;
ids101 = (idsX + 2) + (idsY + (idsZ+1)*GRID.dimY)*GRID.dimX;
ids011 = (idsX + 1) + ((idsY+1)+ (idsZ+1)*GRID.dimY)*GRID.dimX;
ids111 = (idsX + 2) + ((idsY+1)+ (idsZ+1)*GRID.dimY)*GRID.dimX;
%
dess = reshape(dess, [1, GRID.dimX*GRID.dimY*GRID.dimZ]);
tp00 = dess(ids000).*(1-tX) + dess(ids100).*tX;
tp10 = dess(ids010).*(1-tX) + dess(ids110).*tX;
tp01 = dess(ids001).*(1-tX) + dess(ids101).*tX;
tp11 = dess(ids011).*(1-tX) + dess(ids111).*tX;
tp0 = tp00.*(1-tY) + tp10.*tY;
tp1 = tp01.*(1-tY) + tp11.*tY;
vertexDess = tp0.*(1-tZ) + tp1.*tZ;
end