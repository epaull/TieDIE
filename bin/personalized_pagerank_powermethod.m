function x = personalized_pagerank_powermethod(G, start_heats, p)
	if nargin < 3, p = .85; end
	% Eliminate any self-referential links
	G = G - diag(diag(G));
	
	% c = out-degree, r = in-degree
	[n,n] = size(G);
	c = sum(G,1);
	r = sum(G,2);
	% Scale column sums to be 1 (or 0 where there are no out links).
	k = find(c~=0);
	D = sparse(k,k,1./c(k),n,n);

	G = p*G*D;
z = ((1-p)*(c~=0) + (c==0))/n;
%normalize start_heats to 1
start_heats = start_heats/sum(start_heats);
	x = start_heats;
	% here: 
	xprev = 1;
	while sum(abs(xprev-x)) > 0.001
	    xprev = x;
x = G*x + (1-p)*start_heats;
	end
	% Normalize so that sum(x) == 1.
	x = x/sum(x);
