/* 
   This package implements maps which are equivalent
   to combinatorial embeddings of graphs in surfaces. 

   Given a combinatorial graph embedding G -> S, the function
   map_get_presentation computes a topological
   presentation of the fundamental group of the  of the 
*/

/*WARNING: For big n, f==4 will almost never work as 
  most map are disconnected. */
rand_map(n, {f=0})={
	/*differents f's :
		-f=0 : map.
		-f=1 : map with one 
		face.
		-f=2 : reduced map.
		-f=3 : reduced map with
		one face.
		-f=4 : reduced, connected ribbon
		graph.
		-f=5 : reduced, not connected ribbon 
		graph.
		*/
	if(f==0, return([cycs_to_perm(2*n, vector(n, i, [2*i-1,2*i])),rand_perm(2*n)]),
		f==1, return([rand_invol(2*n,n),vectorsmall(2*n,u,(u%(2*n))+1)]),
		f==2, G=rand_map(n,0);
		return(map_reduce(G)),
		f==3, G=rand_map(n,1);
		return(map_reduce(G)),

		f==4, G=map_reduce(rand_map(n,0));
			while(!map_is_connected(G),
				G_=rand_map(n,0);
			);
			return(map_reduce(G)),

		f==5, G=map_connected_components(rand_map(n,0))[1];
			return(G)
	);
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

/*Removes patterns of the form ...aa^-1... in 
faces of G as they are nullhomotopic.*/
/*WARNING: Il y a parfois un bug, s2red est parfois pas
une perm juste un vecsmall, j'ai eu un exemple ou #s2=40
et s2[2]=43 et c'était le seul pb.*/
/*TODO: To be fixed. Currently shouldn't be used.*/
map_reduce(G)={
	if(map_is_reduced(G), return(G));
	my(n,s1,s2);
	[s1,s2]=G;
	n=#s1;

	my(s2inv);
	s2inv=s2^-1;

	my(r,t,t_);
	r=vectorsmall(n,i,i);
	t=vectorsmall(n,i,i);

	/*Patterns of the form ... aa^-1... are
	edges satistying s2(a)=s1(a). First makes
	a and a^-1 into fixed point of s2 and s1.*/
	my(a,b,binv, mark);
	mark=vector(n);
	for(i=1, n,
		if(s1[i]==i, error("s1 has fixed points"));
		if(s2[i]!=s1[i],next);
		if(mark[i], next);
		mark[i]=1;
		mark[s1[i]]=1;

		a=s2inv[i];
		b=i;
		binv=s2[i];

		/*(...a b binv c...)(binv b a)=(...a c...)*/
		t_=maketij(n,i,s1[i]);
		if(a==binv, r=r*t_,
			r=r*make3c(n,binv,b,a)
		);

		t=t*t_;
	);
	my(s1one,s2one);
	s1one=s1*t;
	s2one=s2*r;
	/*
	Then move fixed points to the end and truncate.
	*/
	return(map_reduce(perm_normalize_wrt([s1one,s2one],1)));
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
	to perform a dfs in Gdual and compute a covering tree T.
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

\\Returns a covering tree in G
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
/*a tree in the underlying graph of G^* */
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

		my(g, sone);
		g=(2-(#permcycles(s1*s2)-n/2+#T+1))/2;
		if(g==0,
			sone=vectorsmall(n,i,i);
			s1one=vectorsmall(n,i,i);
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
		
		/*Build s1one and s2one*/
		my(s1one, k, sinv);
		s1one=vectorsmall(n,i,i);
		sone=vectorsmall(n,i,i);
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

/*Act as if there are no punctures at vertices*/
/*of G and lift shortest paths in T from vertex 1*/
/*to each vertex. Then used to build the presentation*/
/*of the fundamental group of the underlying graph of G*/
/*associated to T the covering tree T.*/
map_liftalongT(G, T, TvG)={
	my(s1,s2, s, sinv, n, sc, f);
	[s1,s2]=G;
	s=(s2^-1)*s1;
	sinv=s^-1;
	n=#s1;
	sc=permcycles(s);
	f=#sc;
	vG=map_face_index([s1,s]);


   	my(seedlp, n_one);
	n_one=n-2*(f-1);
	
	my(slpspaths, pointersspaths);
	slpspaths=vector(f);
	pointersspaths=vector(f);
	slpspaths[1]=[-1,0];
	pointersspaths[1]=1;

	my(last, lastvindex, lastvTindex, index, k, e);
	k=1;
	for(i=2, f,
		/*
		   The shortest path from
			the first node of T to the node 
			of index i in the ordering of T.
		*/
		e=T[i-1][1];
		lastvindex=vG[T[i-1][1]];
		lastvTindex=TvG[lastvindex];

		/*Recover index of last instruction in slp*/
		last=pointersspaths[lastvTindex];
		index=n+last;

		slpspaths[k]=[index, e];
		pointersspaths[i]=k;
		k++;
	);

	slpspaths=slpspaths[1..(k-1)];
	return([slpspaths, pointersspaths]);
}

/*Each path gi is homotopic to the unique shortest path in T between
face 1 and i in the dual graph. So that the tree made of the gis 
is isomorphic to T but avoids the hypothetic punctures at the center
of the faces of the dual graph. */
map_liftalongTspecial(G, T, TfG, {type=1})={
	my(s1,s2,s, sinv, n, sc, f);
	[s1,s2]=G;
	if(type, s=s2, s=(s2^-1)*s1);
	sinv=s^-1;
	n=#s1;
	sc=permcycles(s);
	f=#sc;
	fG=map_face_index([s1,s]);

	my(makeslpgammai, maxfacesize = 1);
	/*Encodes a face of Gdual starting at seed (inclusive) and*/
	/*went through in reversed orientation.*/
	makeslpgammai=(u-> 
		my(seed, e, k=1);
		if(u==1,
			seed=1;
		,/*else*/
			seed=T[u-1][2];
		);
		my(slpgammai);
		slpgammai=vector(#sc[fG[seed]]);
		e=seed;
		until(e==seed,
			/*Faut tenir compte des générateurs*/
			if(k==1, 
					slpgammai[k]=[0, s1[e]];
		    ,\\else
					slpgammai[k]=[k-1+n, s1[e]];
				
			);
			e=sinv[e];
			k++;
		);
		maxfacesize=max(maxfacesize, k-1);
		return(slpgammai);
	);
	my(slpsgammai);
	slpsgammai=vector(f,u, makeslpgammai(u));


   	my(seedlp, n_one);
	n_one=n-2*(f-1);
	
	my(slpgis, pointersgi, approxtotlength=maxfacesize*n);
	/*straight line program to compute loopfaces paths.*/
	slpgis=vector(approxtotlength);
	/*Used to read the slp : pointersgi[fTindex]*/
	/*points to the index of the last edge of the path*/
	/*associated to fTindex in slpgis.*/
	pointersgi=vector(f);
	slpgis[1]=[-1,0];
	pointersgi[1]=1;

	my(last, estart, eend, flastindex, flastTindex, k, e);
	k=2;
	for(i=2, f,
		/*
			The path associated to gi is used to join
			the first face to the face 
			of index i in the ordering of T.
		*/
		if(i==2,
			[estart, eend]=[1, T[1][1]];
			last=1;
		,/*else*/
			flastindex=fG[T[i-1][1]];
			flastTindex=TfG[flastindex];
			if(flastTindex==1,
				[estart, eend]=[1, T[i-1][1]];
			,/*else*/
				[estart, eend]=[T[flastTindex-1][2], T[i-1][1]];
			);
			/*Recover index of last instruction in slp*/
			last=pointersgi[flastTindex];
		);

		if(s[estart]==eend,
				/*Add Empty path*/
				/*TODO: transformer en ne rien faire et pas incrémenter k.*/
				slpgis[k]=[0, n+last];
				k++;
		,/*else*/
				e=s[estart];
				until(e==eend,
					/*Unique path contained in face*/
					/*going from s2dual[estart] to */
					/*s2dual^-1[eend] in clockwise */
					/*orientation.*/
					slpgis[k]=[n+last, e];
					e=s[e];
					last=k;
					k++;
				);
		);
		pointersgi[i]=k-1;
	);

	slpgis=slpgis[1..(k-1)];
	return([slpsgammai, slpgis, pointersgi]);
}

/*Assumes G is of genus > 0*/
findab(s2one,s1, seed, index)={
	my(a, ainv, ainvindex);
	a=seed; ainv=s1[a]; ainvindex=index[ainv];
	my(b, binv, binvindex);
	b=s2one[a]; binv=s1[b]; binvindex=index[binv];
	while(binvindex<ainvindex && b!=ainv,
		/*here aindex:=index[a]<bindex:=index[b]<ainvindex by definition*/
		b=s2one[b];
		binv=s1[b];
		binvindex=index[binv];
	);
	return([a,b]);
}

/*
   Performs a single cut and paste on the only cycle of s2one,
 representing the oneface reduction of the graph embedding
 [s1,s2].
*/
cut_and_paste_one(s2one, s1, seedslp, index)={
	my(n, n_one, a,b,vecalpha, vecbeta, vecgamma, vecdelta);
	n=#s1; n_one=permorder(s2one);
	e=seedslp;
	[a,b]=findab(s2one, s1, s2one[e], index);

	my(w);
	w=vector(n_one, u, e=s2one[e]; e);
	e=a;

	my(card);
	card=(index[b]-1)-index[a];
	vecalpha=vector(card, u, e=s2one[e]; e);
	e=s2one[e];
	
	card=(index[s1[a]]-1)-index[b];
	vecbeta=vector(card, u, e=s2one[e]; e);
	e=s2one[e];
	
	card=(index[s1[b]]-1)-index[s1[a]];
	vecgamma=vector(card, u, e=s2one[e]; e);
	e=s2one[e];
	
	vecdelta=vector(n_one-index[s1[b]], u, e=s2one[e]; e);
	
	
	my(slpc, slpd);
	slpd=slpconcat([vectoslp(vecalpha,n),[[0,b]], vectoslp(vecbeta,n)],n,,1);
	slpc=slpconcat([slpinvert(vectoslp(vecbeta,n),s1,n), slpinvert(vectoslp(vecgamma,n),s1,n), [[0, a]] ],n,,1);

	\\updated_w=gamma beta c d c^-1d^-1 alpha delta
	\\ a and b are used because they are not in the vecs
	\\ the side pairing then associates to a and b
	\\ the evaluation of slpc and slpd respectively
	my(updated_w, updated_s2one, updated_index);
	updated_w=concat([vecgamma, vecbeta, [a], [b], [s1[a]], [s1[b]], vecalpha, vecdelta]);
	updated_s2one=cycs_to_perm(n, [updated_w]);
	updated_index=map_indexcycle(updated_s2one, updated_w[1]);
		
	return([a,b, vecalpha, vecbeta, vecgamma, vecdelta, slpc, slpd, w, updated_w, updated_s2one, updated_index]);
}

/* 
   Used to filter 2g generators and build a relation in buildpres
   function in map.gp.
 */
buildrel_and_pointers(pointers, seed, s, s1, index, seen)={
	if(#pointers==1, return([pointers,[]]));
	my(e=seed, rel, j=1, newpointers, newindex, m=0);
	newindex=index;
	newpointers=vector(#pointers/2);
	rel=vector(#pointers);

	for(i=1, permorder(s),
		if(seen[s1[e]], 
			rel[i]=-newindex[s1[e]];
			m++;
		,/*else*/
			newpointers[j]=pointers[i];
			j++;
			newindex[e]+=-m;
			rel[i]=newindex[e];
		);
		e=s[e];
	);
	\\error("");
	return([newpointers, rel]);
}

map_surface_presentation(Gone, seedslp, {type="oneword"})={
	my(s1one, s2one,s2oneinv);
	[s1one,s2one]=Gone;
	n=#s1one;
	s2oneinv=s2one^-1;
	n_one=permorder(s2one);
	if(n_one==1, return([[],[],[]]));
	f=(n-n_one)/2+1;

	my(slp, pointers, seen, index, e, w);
	seen=vectorsmall(n);
	w=vectorsmall(n_one);
	e=seedslp;
	index=map_indexcycle(s2one,s2one[e]);
	\\Each edge e points to e(1).
	for(i=1, n_one,
		e=s2one[e];\
		if(!seen[s1one[e]], seen[e]=1);\
		w[i]=e;
	);
	e=seedslp;

	my(vecalpha, vecbeta, vecgamma, vecdelta,\
		slpalpha, slpbeta, slpgamma, slpdelta,\
	   	   	slpc,slpd);
	my(a,b, updated_w, updated_s2one, updated_index);
	if(type=="oneword",
		pointers=vectorsmall(n_one, i, i);
		slp=vector(n_one, i,\
			   	[0, w[i]]);
	,type=="onehandle",/*else if*/
		\\w_one=gamma beta c d c^-1d^-1 alpha delta
		\\ Recover updated pointers of gis before concatening
		\\ the loopfaces to slp and pointers.
		\\ In building pointers we wish to point to the 2*g
		\\ generators and the f loopfaces. As each edge
		\\ appears with its inverse in the one word, we
		\\ point only to the first of the two appearing
		\\ following the w_one.
		[a,b, vecalpha, vecbeta,vecgamma, vecdelta, slpc, slpd,
	   	w, updated_w, updated_s2one, updated_index]=\
			cut_and_paste_one(s2one,s1one,seedslp,index);

		slpalpha=vectoslp(vecalpha, n, 0);
		slpbeta=vectoslp(vecbeta, n, 0);
		slpgamma=vectoslp(vecgamma, n, 0);
		slpdelta=vectoslp(vecdelta, n, 0);

		[slp, pointers]=slpconcat(
			\\ SLPS
			[slpgamma, slpbeta, slpc, slpd], n,\
			\\ POINTERS
			[vector(#slpgamma, i,
				if(!seen[s1one[slpgamma[i][2]]], seen[slpgamma[i][2]]=1); i),\
			vector(#slpbeta, i,
				if(!seen[s1one[slpbeta[i][2]]], seen[slpbeta[i][2]]=1); i),\
			[#slpc], [#slpd]]\
		);
		
		\\ Add slpc^-1 and slpd^-1
		slp=concat([slp,\
			[[-(n+#slpgamma+#slpbeta+#slpc), -1]],\
		   	[[-(n+#slpgamma+#slpbeta+#slpc+#slpd), -1]]]);
		pointers=concat([pointers,[n+#slpgamma+#slpbeta+#slpc+#slpd+1,n+#slpgamma+#slpbeta+#slpc+#slpd+2]]);
		
		[slp, pointers]=slpconcat(
			\\ SLPS
			[slp,\
	   		slpalpha, slpdelta], n,\
			\\ POINTERS
			[pointers,
			vector(#slpalpha, i,
				if(!seen[s1one[slpalpha[i][2]]], seen[slpalpha[i][2]]=1); i),\
			vector(#slpdelta, i,
				if(!seen[s1one[slpdelta[i][2]]], seen[slpdelta[i][2]]=1); i)]\
		);
		s2one=updated_s2one;
		index=updated_index;
		\\ a and b are now c and d respectively for updated_index
		\\ This works as only their index in updated_w is used.
		seen[a]=1;
		seen[b]=1;
	);
	\\Pointers has size 4*g at this point, that is it points on
	\\2*g topological generators and their inverses. 
	\\This filters 2*g topological generators and builds
	\\the appropriate relation.
	my(rel);
	[pointers,rel]=buildrel_and_pointers(pointers, slp[1][2], s2one, s1one, index, seen);
	return([slp, pointers, rel]);
}

map_loopfaces(n, slpsgammai, slpgis, pointersgi)={
	my(f, n_one, rel);
	f=#slpsgammai;
	n_one=n-2*(f-1);
	
	rel=vector(f, i, n_one/2+i);
	\\ Makes slps of gammais into a single slp
   	my(slpgammais, pointersgammai);
	[slpgammais, pointersgammai]=slpconcat(slpsgammai, n,\
		vector(#slpsgammai, u, [#slpsgammai[u]])\
	);
	my(slp, pointers);
	\\ Add slp of gis and gammais 
	[slp, pointers]=slpconcat(
		\\ SLPS
		[slpgis, slpgammais], n,\
		\\ POINTERS
		[pointersgi, pointersgammai]\
	);

	\\ Add loopfaces, slp with generator set <gi, gammaj;i,j>
	\\ of size #pointers==2*f
	my(slp_loopfaces, pointers_loopfaces);
	slp_loopfaces=concat(\
			\\SLPS
			vector(f, u,\
			\\ u-th gammai
			[[u, f+u],\
			\\ compute gi^-1
			[-u, -1],\
			\\ (gi.gammai).(gi^-1)
			[(2*f+3*u)-2, (2*f+3*u)-1]])\
		);
	pointers_loopfaces=vector(f,u, 3*u);
	[slp,pointers]=slpcompose([n, 2*f], [slp, slp_loopfaces],\
			[pointers, pointers_loopfaces]
	);
	return([slp, pointers, rel]);
}

map_topological_presentation_fromT(G, T, TfG, type)={
	my(Gone, s1one, s2one, seed);
	[s1one, s2one, seed]=map_gluealongT(G, T);
	Gone=[s1one,s2one];
	my(slpsgammai, slpgis, pointersgi);
	[slpsgammai, slpgis, pointersgi]=map_liftalongTspecial(G,T, TfG);

	my(surface_slp, surface_pointers,surface_rel,\
			loopfaces_slp, loopfaces_pointers, face_rel);
	[surface_slp, surface_pointers, surface_rel]=\
			map_surface_presentation(Gone, seed,  type);
	[loopfaces_slp, loopfaces_pointers, face_rel]=\
			map_loopfaces(#G[1], slpsgammai, slpgis, pointersgi);

	my(fullslp, pointers, rels);
	\\ Concatenate without connecting
	[fullslp,pointers]=slpconcat(\
			\\ SLPS
			[surface_slp, loopfaces_slp], #G[1],\
			\\ POINTERS
			[surface_pointers, loopfaces_pointers]);
	rel=concat(surface_rel, face_rel);
	return([fullslp, pointers, rel]);
}

/*Recover presentation of a map G -> S*/
/*associated to its dual map G^*->S*/
map_topological_presentation(G, {type="oneword"})={
	my(Gdual);
	Gdual=map_dual(G);

	my(T, TfGdual);
	[T, TfGdual]=map_getT(Gdual);
	my(slp, pointers, rel);
	[slp, pointers, rel]=map_topological_presentation_fromT(Gdual, T, TfGdual, type);
	return([slp, pointers, rel, T, TfGdual]);
}

/*Returns a side pairing for the map G/T->S as described in framework.txt*/
/*Here G should be tought as G^*->S for a map G-> S*/
map_quotientT(G, T, TvG)={
	my(vG);
	vG=map_vertex_index(G);

	my(s0one, s1one, s2one, seed);
	[s1one, s2one, seed]=map_gluealongT(G, T, 0);
	s0one=(s2one^-1)*s1one;

	my(slpspaths, pointersspaths);
	/*Slp and pointers of shortest paths*/
	[slpspaths, pointersspaths]=map_liftalongT(G, T, TvG);

	my(slpold, pointersold);
	/*Slp and pointers of old side pairing*/
	[slpold, pointersold]=map_surface_presentation([s1one, s0one], seed, "oneword");

	my(slp_conjs, pointers_conjs, e);
	e=seed;
	/*Slp and pointers of new side pairing*/
	slp_conjs=concat(\
			\\SLPS
			vector(n_one, u,\
			\\ c.e
			e=s0one[e];
			[[u, f+TvG[vG[e]]],\
			\\ compute c^-1
			[-u, -1],\
			\\ (c.e).(c-1)
			[(n_one+f+3*u)-2, (n_one+f+3*u)-1]])\
		);
	pointers_conjs=vector(n_one,u, 3*u);
	my(slp, pointers);
	[slp, pointers]=slpconcat(
			\\SLPS
			[slpold, slpspaths], n,\
			\\POINTERS
			[pointersold, pointerspaths]\
		);

	/*Slp and pointers of new side pairing*/
	[slp, pointers]=slpcompose([n,n_one+f], [slp, slp_conjs],[pointers, pointers_conjs]);
	return([slp, pointers, s1one, s2one, seed]);
}
map_quotient_spair(G)={
	my(T, TvG);
	[T, TvG]=map_getT(G, 0);
	return(map_quotientT(G, T, TvG));
}
\\ Build a covering map coming from a monodromy
\\ action.
\\ The new edge set is Erev=E x {1,..., d}, we view it as 
\\ {1,...,#E*d} with lexicographic ordering.
\\ The covering vertex permutation s0rev is s0 x id_{1,...,d}
\\ while s1rev is s1(_) x (monodromy(_)(s1(_))).
\\ TODO: Le sidepairing devrait être calculé dans fdompres.
map_from_monodromy(G, d, monodromy)={
	my(s0, s1, s2, n);
	[s1,s2]=G;
	s0=(s2^-1)*s1;
	n=#s0;

	my(s0rev, s1rev, m=d*n);
	s0rev=vectorsmall(m);
	s1rev=vectorsmall(m);
	for(j=1, d,
		for(i=1, n,
			s0rev[i+(j-1)*n]=s0[i]+(j-1)*n;
			s1rev[i+(j-1)*n]=s1[i]+(monodromy[s1[i]][j]-1)*n;
		);
	);
	return([s1rev, s1rev*s0rev^-1]);
}

\\ Given a map G->S with orbifold points at vertices of G
\\ and a covering orbifold Srev -> S coming from monodromy.
\\ Computes a one face map Grevone -> Srev with orbifold
\\ points at its vertices as well as a side pairing.
map_monodromy_getonefacespair(G, d, monodromy)={
	my(Gdual, Gdualrev, Gdualrevdual);
	Gdual=map_dual(G);
	Gdualrev=map_from_monodromy(Gdual, d, monodromy);


	my(s0one, s1one, s2one, index);
	\\ To obtain the side pairing associated to s1one, order
	\\ the edges in the support (the only cycle) of s0one starting
	\\ at seed in clock-wise orientation (s0one). Then the i'th
	\\ edge, that is e s.t index[e]=i is associated to slp[pointers[i]].
	[slp, pointers, s1one, s2one, seed]=map_quotient_spair(Gdualrev);
	index=map_indexcycle(s0one, s0one[seed]);

	my(s0,s1,s2,g);
	\\ Here s1=g^-1.s1one.g and s2=g^-1.s2one.g so that
	\\ the side pairing of [s1,s2] is obtained from
	\\ that of [s1one, s2one] by letting e in s1 associate
	\\ slp[pointers[i]] if i=index[g[e]]
	[s1, s2, g]=perm_normalize_wrt([s1one,s2one],1,,1);
	g=g^-1[1..#s1];

	normalizer=vector(#s1, e, index[g[e]]);
	slp=slpnormalize(slp, normalizer, #s1one);

	my(Gdualrevone, Grevone);
	Gdualrevone=[s1, s2];
	s0=(s2^-1)*s1;
	Grevone=[s1, s0];

	return([Grevone, slp, pointers]);
}
