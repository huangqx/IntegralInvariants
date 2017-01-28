function [] = save_descriptors(vertexDess, filepath)
%
vertexDess = double(vertexDess');
save(filepath, 'vertexDess', '-ascii');