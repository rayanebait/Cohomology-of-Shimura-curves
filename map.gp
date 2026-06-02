/* 
   This package implements maps which are equivalent
   to combinatorial embeddings of graphs in surfaces. 
*/

/*WARNING: For big n, most maps are disconnected. */
rand_map(n, {f=0})={
	/*differents f's :
		-f=0 : map.
		-f=1 : map with one 
		face.
		*/
	if(f==0, return([cycs_to_perm(2*n, vector(n, i, [2*i-1,2*i])),rand_perm(2*n)]),
		f==1, return([rand_invol(2*n,n),vectorsmall(2*n,u,(u%(2*n))+1)]),
	);
	return();
}

/*Given an involution representing a side pairing of a fundamental
  domain, with possible fixed points, slices each fixed edge in 2
  and returns a map together with a map associating to
  each edge of the map the index of its associated 
  generator.*/
map_from_invol(invol)={
	my(n,sliced,k,ki);
	n=#invol;
	k=0;
	/*i becomes i+ki[i]*/
	ki=vector(n);
	sliced=vector(n);
	for(i=1,n,
		ki[i]=k;
		if(invol[i]!=i, next);
		k++;
		sliced[i]=1;
	);

	/* If e is a fixed point of invol, slice e into
	a and a^-1, geometrically add a vertex in 
	between e(0) and e(1). Then slide the indexes
	on the right */
	/*NOTE : l'arête e correspond à g_e dans elts*/
	/*à la fin on veut évaluer des mots d'arêtes */
	/*associés */
	
	my(s1,s2);
	s1=vectorsmall(n+k,i,i);
	s2=perm_iplusk(n+k,1);

	my(inew,iinvnew, g);
	g=vectorsmall(n+k);
	for(i=1,n,
		inew =i+ki[i];
		g[inew]=i;
		if(sliced[i],
			s1[inew]=inew+1;
			s1[inew+1]=inew;
			g[inew+1]=i;
			next;
		);
		iinvnew=invol[i]+ki[invol[i]];
		s1[inew]=iinvnew;
		s1[iinvnew]=inew;
	);

	return([[s1,s2],g]);
}


map_dual(G)={
	my(s0,s1,s2, Gdual);
	[s1,s2]=G;
	s0=s2^-1*s1;
	Gdual=[s1, s1*s0*s1];
	return(Gdual);
}

/*Return the genus, number of edges vertices and faces of the
  graph embedding G*/
map_numbers(G)={
	if(#G[2]<=1, return([0,0,0,0]));
	my(e,f,v,g);
	e=#permcycles(G[1]);
	f=#permcycles(G[2]);
	v=#permcycles(G[2]^(-1)*G[1]);

	g=(v-e+f-2)\(-2);
	return([v,e,f,g]);
}

/*normalizes s2 so that 

  s2=(1...f1)(f1+1...f1+f2)...((f1+...+f_f-1) ... f1+...+f_f)

where f=#permcycles(s2). */
map_normalize(G)={
	my(n, s1,s2);	
	[s1,s2]=G;
	n=#s1;
	
	my(s2c,g, f, k, fi, cardfi);
	s2c=permcycles(s2);
	g=vectorsmall(n,u,u);

	f=#s2c;
	k=1;
	for(i=1,f,
		fi=s2c[i];
		cardfi=#fi;
		for(j=1, cardfi,
			g[fi[j]]=k;
			k=k+1;
		);
	);
	my(s1one,s2one);
	s2one=permconj(s2,g);
	s1one=permconj(s1,g);
	return([s1one,s2one]);
}

/*Checks for any patterns of the form ...aa^-1...
in the faces of G.*/
map_is_reduced(G)={
	my(s1,s2,n,a);
	[s1,s2]=G;
	n=#s2;
	if(n<=1, return(1));
	a=1;

	for(i=1, n,
		a=a*(s1[i]!=s2[i]);
	);
	return(a);
}

/*Given the ordering induced by permcycles on cycles of s2 (faces),
associates to each oriented edge the index of the face it belongs to.*/
map_face_index(G)={
	my(s2, n, fG);
	s2=G[2];
	n=#s2;
	fG=vector(n);

	my(s2c,face);
	s2c=permcycles(s2);
	for(f=1,#s2c, 
		face=s2c[f];
		/*To each face associate an index i,
		to each edge of the face associate
		the this index i.*/
		foreach(Vec(face), i,
			fG[i]=f;
		);
	);
	return(fG);
}


/*
Equivalent to map_face_index([s1,s0]).
*/
map_vertex_index(G)={
	my(s0, s1, s2, n, vG, vertices);

	[s1,s2]=G;
	s0=s2^-1*s1;
	n=#s2;
	vG=vectorsmall(n);

	vertices=permcycles(s0);
	for(v=1,#vertices, 
		vertex=vertices[v];
		foreach(Vec(vertex), i,
			vG[i]=v;
		);
	);
	return(vG);
}


/*NOTE: The set of edges of a map G on edges 
E={1,...,n} is given by :
	-E(G)=E
	-V(G)={(ci, i) | i=1,...,vG} where 
	Prod_i c_i = s0 is the cycle decomposition 
	of s0 and the indexing is the one computed in  
	map_vertex_index. Recall that s0 is a 
	permutation of E so that it makes sense to ask
	if e is in a given cycle ci.
  This function returns a map Phi : E(G)->V(G)xV(G) such
  that Phi(e)=(e(0),e(1)):=(e,s1[e]).*/
/*NOTE:O(n^2) space as each cycle is O(n) and |E|=n.
Here only for completeness.*/
map_graph(G)={
	my(vG, s1, s2, Phi, n);
	vG=map_vertex_index(G);
	[s1, s2]=G;
	n=#s1;


	Phi = vector(n);
	for(i=1, n,
		Phi[i]=[vG[i],vG[s1[i]]];
	);
	return(Phi);
}

/*NOTE: Standard depth first search using the underlying
graph
of G output by map_graph. Currently O(n^2) time and
space due to the representation of G as a graph. Can
be done in much less using a vertex index and s0.*/
map_dfsgen(G)={
	my(s0, cardvertices);
	s0=G[2]^-1*G[1];
	cardvertices=#permcycles(s0);

	my(PhiG, v, explored, to_explore, e);
	PhiG=map_graph(G);
	v=PhiG[1][1];
	explored=vector(cardvertices);
	to_explore=List();

	listput(~to_explore,v);
	/*Counts explored*/
	e=0;
	while(#to_explore!=0, 
		c=#to_explore;
		v=to_explore[c];
		listpop(~to_explore);
		explored[v[2]]=1;
		e=e+1;
		foreach(Vec(v[1]), j,
			[vj,vjinv]=PhiG[j];
			if(!explored[vjinv[2]], listput(~to_explore, vjinv););
		);
	);
	return([explored,e]);
}


map_dfsinit(G, ~Phi, ~data, {type=1})={
	my(s,vG,s1,sc);
	s1=G[1];
	if(type, s=G[2], s=(G[2]^-1)*s1);
	vG=map_face_index([s1,s]);
	/*Go through the faces in reversed orientation. Needed for the slps.*/
	s=s^-1;
	sc=permcycles(s);

	my(explored, vdfsindex, dfsvG, seed);
	explored=vector(#sc);
	vdfsindex=1;
	dfsvG=vector(#sc);
	seed=1;

	my(bufPhi);
	bufPhi=[s1, s, sc, vG];
	for(i=1, #Phi, Phi[i]=bufPhi[i]);

	my(bufdata, tempT);
	tempT=List();
	bufdata=[explored, tempT, dfsvG, vdfsindex, seed, vG];
	for(i=1, #data, data[i]=bufdata[i]);

	return();
}

/*
	Performs a depth first search in the graph associated to
	a graph embedding G starting at edge seed=1.

	Initialize my(data) as a reference with map_dfs(G, ~data)
	to perform a dfs in Gdual and compute a spanning tree T.
	Phi has the form
		[s1, s, sc, vG]=Phi;
	data has the form
		[explored, T, dfsvG, vdfsindex, seed]=data
 */
/*type=0 does a dfs in the vertices of G while type=1 in the*/
/*faces of G, equivalently in the vertices of the dual graph.*/
map_dfs(~Phi, ~data, {type=1})={
	my(G);
	/*rec_prof++;*/
	/*if(rec_prof > 5, breakpoint());*/
	if(#Phi==2,
		G=Phi;
		Phi=vector(4);
		map_dfsinit(G, ~Phi, ~data, type);
	,/*else*/
		/*data[4] is the index of the current recursion*/
		data[4]++;
	);

	/*Data related to current vertex.*/
	my(seed, v, vindex);
	/*seed is the first edge to be visited in this*/
	/*vertex*/
	seed=data[5];
	vindex=Phi[4][seed];
	v=Phi[3][vindex];


	/*Map from the permcycles ordering to*/
	/*the dfs ordering of the vertices. Used*/
	/*mainly with the straight line program (slp)*/
	data[3][vindex]=data[4];
	/*data[3] associates to a vertex index its */
	/*index for the dfs ordering.*/


	/*Mark current vertex*/
	data[1][vindex]=1;
	if(data[4]==1 && #v==1, 
		/*Can happen only at first iteration */
		listput(~data[2],[seed,Phi[1][seed]]);
		/*Update seed*/
		data[5]=Phi[1][seed];

		map_dfs(~Phi,~data, type);
		return();
	);
	
	my(j, jinv, vjinvindex);
	j=seed;
	until(j==seed,
		jinv=Phi[1][j];
		vjinvindex=Phi[4][jinv];


		/*Already explored*/
		if(data[1][vjinvindex],
			j=Phi[2][j];
			next;
		);
		/*Unexplored*/

		listput(~data[2],[j,Phi[1][j]]);
		/*Update seed*/
		data[5]=jinv;
		map_dfs(~Phi, ~data, type);

		/*Increment*/
		j=Phi[2][j];
	);
	/*rec_prof--;*/
	return();
}

\\Returns a spanning tree in G
\\in the format of a vector of couple
\\ of edges [e, s1(e)] such that 
\\e and s1(e) lie in different vertices of G
\\It also returns the ordering of
\\ permcycles(s2) coming from the dfs.
\\ type is as in map_dfs
map_getT(G, {type=1})={
	my(data=vector(6), T, TfG);
	map_dfs(~G, ~data,type);
	[T, TfG]=data[2..3];
	return([Vec(T),TfG]);
}

map_is_connected(G)={
	if(#permcycles(G[2])<=1, return(1));
	my(data=vector(6), explored);
	map_dfs(~G,~data);
	explored=data[1];
	foreach(explored,i, if(!i, return(0)));
	return(1);
}

/*Utility function : Builds a list of maps CC such that
  each is isomorphic to a connected component of G.*/
map_connected_components(G, {CC=List()})={
	my(n, s1, s2, s2c, explored, e);
	if(#G[1]==0 || #G[2]==0, return(CC));
	my(data=vector(6));
	map_dfs(~G,~data);
	explored = data[1];
	e=0;
	foreach(explored, i, e=e+i);

	[s1,s2]=G;
	n=#s1;
	s2c=permcycles(s2);

	/*NOTE: Ending condition*/
	if(e==#s2c, listput(~CC, G);
		return(CC);
	);

	/*NOTE: Build a connected component 
	aswell as a map with the remaining 
	connected components.*/
	my(s2CCc, s2leftc, s1CC, s1left);
	s2CCc=List();
	s2leftc=List();
	s1CC=vectorsmall(n,i,i);
	s1left=vectorsmall(n,i,i);

	my(c);
	for(i=1, #explored,
		c=Vec(s2c[i]);
		if(explored[i],
			/*Build map out of 
			the connected component*/
			listput(~s2CCc, c); 
			foreach(c, j,
				s1CC[j]=s1[j];
			);
		,/*else*/
			/*Build map out of 
			the remaining
			connected components*/
			listput(~s2leftc, c);
			foreach(c, j,
				s1left[j]=s1[j];
			);
		);
	);
	my(s2CC, s2left);
	s2CC=cycs_to_perm(n, s2CCc);
	s2left=cycs_to_perm(n, s2leftc);

	listput(~CC, perm_normalize_wrt([s1CC,s2CC],1));

	my(Gleft);
	Gleft=perm_normalize_wrt([s1left,s2left],1);

	/*NOTE: Recursive call on the map left*/
	CC=map_connected_components(Gleft, CC);

	return(CC);
}

/*Builds the map obtained by gluing the faces of G along T*/
/*a tree in the underlying graph of the dual graph G^*. */
map_gluealongT(G, T, {type=1}, {withGone=0})={
		/*Return if it already has one face that is T is empty*/
		if(#T==0, 
			if(withGone,
				if(type, 
					return([G[1],G[2], 1, G]);
				,/*else*/
					return([G[1], (G[2]^-1)*G[1], 1, G]);
				);
			,/*else*/
				if(type,
					return([G[1], G[2], 1]);
				,/*else*/
					return([G[1], (G[2]^-1)*G[1], 1]);
				);
			);
		);

		my(n,s1,s2);
		[s1,s2]=G;
		n=#s1;

		my(g, s1one, sone);
		g=(2-(#permcycles(s1*s2)-n/2+#T+1))/2;
		s1one=vectorsmall(n,i,i);
		sone=vectorsmall(n,i,i);
		if(g==0,
			if(withGone,
				return([s1one,sone,1, [[],[]]]);
			,/*else*/
				return([s1one,sone,1]);
			);
		);
		
		my(in_T);
		in_T=vector(n);
		
		foreach(T, e,
			in_T[e[1]]=1;
			in_T[e[2]]=1;
		); 
		
		my(sinv);
		/*Build s1one and s2one*/
		if(type,
			sinv=s2^-1;
		,/*else*/
			sinv=s1*s2;
		);
		for(i=1, n, 
			if(in_T[i], next);
			s1one[i]=s1[i];
			k=sinv[i];
			\\k=s2[i];
			while(in_T[k],
				\\k=s2[s1[k]];
				k=sinv[s1[k]];
			);
			sone[i]=k;
		);
		sone=sone^-1;

		my(seed=1);
		while(in_T[seed],
			seed=sinv[s1[seed]];
		);
		if(!type, sone=s1one*sone^-1);
		if(withGone, 
			return([s1one, sone, seed, perm_normalize_wrt([s1one,sone], 1, in_T)]);
		,/*else*/
			return([s1one, sone, seed]);
		);
}

