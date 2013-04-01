
disp('Loading the network from file /tmp/in.tab')
fid = fopen('/tmp/in.tab');
S = textscan(fid,'%s %s', 'Delimiter', '\t');
fclose(fid);

disp('Loading the heat sources vector from file /tmp/heats.tab')
fidh = fopen('/tmp/heats.tab');
heat_sources = textscan(fidh,'%s %f', 'Delimiter', '\t');
fclose(fidh);


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
diagIdx      = sub2ind(size(Adj), 1:length(Adj), 1:length(Adj));
Adj(idx1)    = 1;
Adj(diagIdx) = 0;

% fix input heat vector
keys = Nodes;
heat_vec = zeros(length(keys), 1);
for j=1:length(keys)
   % Indices in the user-supplied vector.
   k = find(strcmp(keys(j),heat_sources{1}));

   if(length(k) > 1)
      disp('ERROR: Repeat identifiers in heat sources file; 1st column must have unique ids.')
      return;
   elseif(length(k) == 1)
      heat_vec(j) = heat_sources{2}(k);
   end
end

% use 0.85 as a default. 
disp('Writing result heat vector to /tmp/result.tab');
result = personalized_pagerank_powermethod(Adj, heat_vec, 0.85);
result_str = strcat(Nodes, ':', num2str(result));
fid=fopen('/tmp/result.tab','wt');
for i=1:length(result_str)
	fprintf(fid,'%s\n',result_str{i});
end
fclose(fid);
exit;
