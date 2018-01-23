

% Time of evolution of the network. This is taken from the HotNet
% paper, Vandin et al. PSB (2012). Qi et al (Bader lab) in Genome
% Research (2008) report values of 0.1 to 2 gave equally good
% prediction using PPI nets. They tried values between 0.001 up
% to 256.
t = %TIME%;

fid = fopen('/tmp/in.tab');
S = textscan(fid,'%s %s', 'Delimiter', '\t');
fclose(fid);

% Convert the network edges into indexes.
Links = [S{1},S{2}];
[Nodes,dummy,idx] = unique(Links);

disp('Converting node identifiers to indices.')
E = zeros(size(Links));
for j=1:length(Nodes)
   E = E + j*strcmp(Nodes(j), Links);
end

% Adj is the adjacenceny matrix
disp('Computing the adjacency matrix.')
Adj          = sparse(1:length(Nodes), 1:length(Nodes), 0 );
idx1         = sub2ind(size(Adj), E(:,1),E(:,2));
idx2         = sub2ind(size(Adj), E(:,2),E(:,1));
diagIdx      = sub2ind(size(Adj), 1:length(Adj), 1:length(Adj));
Adj(idx1)    = 1;
Adj(idx2)    = 1;
Adj(diagIdx) = 0;

% L is the graph Laplacian
disp('Computing the graph Laplacian.')
L = diag(sum(Adj,1)) - Adj;

% K is the heat diffusion kernel of the graph.
K = [];
disp(sprintf('Calculating the heat diffusion kernel at time t=%f', t))
K = expm(-t*L);
disp(sprintf('Saving the heat diffusion kernel to file %OUTPUT%'))

% old MATLAB 2012 version
% printmatrix(K, char(Nodes), char(Nodes), '%OUTPUT%', 'w');

% Open a file for writing
fid = fopen('%OUTPUT%','w');

% Write the first row: node names
fprintf(fid,'Key')
for ii = 1:size(Nodes,1)
    fprintf(fid,'\t%s',Nodes{ii});
end
fprintf(fid,'\n');

% Write the heat matrix
for ii = 1:size(K,1)
    fprintf(fid,'%s',Nodes{ii});
    fprintf(fid,'\t%g',K(ii,:));
    fprintf(fid,'\n');
end

fclose(fid);
exit;

% disp(sprintf('Saving the heat diffusion kernel to binary file kernel.mat'))
% save 'kernel.mat' K Nodes
