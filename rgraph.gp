/* 
   This package implements ribbon graphs which are equivalent
   to combinatorial embeddings of graphs in surfaces. 

   Given a combinatorial graph embedding G -> S, the function
   rgraph_get_presentation computes a topological
   presentation of the fundamental group of the  of the 
*/

/*WARNING: For big n, f==4 will almost never work as 
  most ribbon graph are disconnected. */
rand_rgraph(n, {f=0})={
	/*differents f's :
		-f=0 : ribbon graph.
		-f=1 : ribbon graph with one 
		face.
		-f=2 : reduced ribbon graph.
		-f=3 : reduced ribbon graph with
		one face.
		-f=4 : reduced, connected ribbon
		graph.
		-f=5 : reduced, not connected ribbon 
		graph.
		*/
	if(f==0, return([cycles_to_perm(2*n, vector(n, i, [2*i-1,2*i])),rand_perm(2*n)]),
		f==1, return([rand_invol(2*n,n),vectorsmall(2*n,u,(u%(2*n))+1)]),
		f==2, G=rand_rgraph(n,0);
		return(rgraph_reduce(G)),
		f==3, G=rand_rgraph(n,1);
		return(rgraph_reduce(G)),

		f==4, G=rgraph_reduce(rand_rgraph(n,0));
			while(!rgraph_is_connected(G),
				G_=rand_rgraph(n,0);
			);
			return(rgraph_reduce(G)),

		f==5, G=rgraph_connected_components(rand_rgraph(n,0))[1];
			return(G)
	);
}

/*Given an involution representing a side pairing of a fundamental
  domain, with possible fixed points, slices each fixed edge in 2
  and returns a ribbon graph together with a map associating to
  each edge of the ribbon graph the index of its associated 
  generator.*/
rgraph_from_invol(invol)={
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


rgraph_dual(G)={
	my(s0,s1,s2, Gdual);
	[s1,s2]=G;
	s0=s2^-1*s1;
	Gdual=[s1, s1*s0*s1];
	return(Gdual);
}

/*Return the genus, number of edges vertices and faces of the
  graph embedding G*/
rgraph_numbers(G)={
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
rgraph_normalize(G)={
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
	my(s1_,s2_);
	s2_=permconj(s2,g);
	s1_=permconj(s1,g);
	return([s1_,s2_]);
}

/*Removes patterns of the form ...aa^-1... in 
faces of G as they are nullhomotopic.*/
/*WARNING: Il y a parfois un bug, s2red est parfois pas
une perm juste un vecsmall, j'ai eu un exemple ou #s2=40
et s2[2]=43 et c'était le seul pb.*/
/*TODO: To be fixed. Currently shouldn't be used.*/
rgraph_reduce(G)={
	if(rgraph_is_reduced(G), return(G));
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
	my(s1_,s2_);
	s1_=s1*t;
	s2_=s2*r;
	/*
	Then move fixed points to the end and truncate.
	*/
	return(rgraph_reduce(perm_normalize_wrt([s1_,s2_],1)));
}

/*Checks for any patterns of the form ...aa^-1...
in the faces of G.*/
rgraph_is_reduced(G)={
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
rgraph_face_index(G)={
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
Equivalent to rgraph_face_index([s1,s0]).
*/
rgraph_vertex_index(G)={
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


/*NOTE: The set of edges of a ribbon graph G on edges 
E={1,...,n} is given by :
	-E(G)=E
	-V(G)={(ci, i) | i=1,...,vG} where 
	Prod_i c_i = s0 is the cycle decomposition 
	of s0 and the indexing is the one computed in  
	rgraph_vertex_index. Recall that s0 is a 
	permutation of E so that it makes sense to ask
	if e is in a given cycle ci.
  This function returns a map Phi : E(G)->V(G)xV(G) such
  that Phi(e)=(e(0),e(1)):=(e,s1[e]).*/
/*NOTE:O(n^2) space as each cycle is O(n) and |E|=n.
Here only for completeness.*/
rgraph_graph(G)={
	my(vG, s1, s2, Phi, n);
	vG=rgraph_vertex_index(G);
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
of G output by rgraph_graph. Currently O(n^2) time and
space due to the representation of G as a graph. Can
be done in much less using a vertex index and s0.*/
rgraph_dfsgen(G)={
	my(s0, cardvertices);
	s0=G[2]^-1*G[1];
	cardvertices=#permcycles(s0);

	my(PhiG, v, explored, to_explore, e);
	PhiG=rgraph_graph(G);
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


/*Performs a depth first search in the graph associated to
a graph embedding G starting at edge seed=1.

Initialize my(data) as a reference with rgraph_dfs(G, ~data)
to perform a dfs in Gdual and compute a covering tree T.
Phi has the form
		[s1, s, sc, vG]=Phi;
data has the form
		[explored, T, slp, dfsvG, vdfsindex, seed]=data


At the end data has the form : [explored,T,dfsvG,fdfsindex,seed, vG]=data*/
rgraph_dfs(~Phi, ~data)={
	/*Data of G*/
	my(s1, s, sc, vG, T, seed);
	/*rec_prof++;*/
	/*if(rec_prof > 5, breakpoint());*/
	if(#Phi==2,
		/*Initialize*/
		my(G);
		G=Phi;
		/*Dfs is on Gdual.*/
		s=G[2];
		vG=rgraph_face_index(G);

		/*Go through the faces in reversed orientation. Needed for the slps.*/
		s=s^-1;
		s1=G[1];
		sc=permcycles(s);

		my(explored, vdfsindex, dfsvG);
		explored=vector(#sc);
		vdfsindex=1;
		dfsvG=vector(#sc);
		seed=1;

		my(tempT);
		tempT=List();

		Phi=[s1, s, sc, vG];
		my(bufdata);
		bufdata=[explored, tempT, dfsvG, vdfsindex, seed, vG];
		for(i=1, #data, data[i]=bufdata[i]);
	,/*else*/
		/*data[4] is the index of the current recursion*/
		data[4]++;
	);

	/*Data related to current vertex.*/
	my(v, vindex);
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

		rgraph_dfs(~Phi,~data);
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
		rgraph_dfs(~Phi, ~data);

		/*Increment*/
		j=Phi[2][j];
	);
	/*rec_prof--;*/
	return();
}

rgraph_is_connected(G)={
	if(#permcycles(G[2])<=1, return(1));
	my(data=vector(6), explored);
	rgraph_dfs(~G,~data);
	explored=data[1];
	foreach(explored,i, if(!i, return(0)));
	return(1);
}

/*Utility function : Builds a list of ribbon graphs CC such that
  each is isomorphic to a connected component of G.*/
rgraph_connected_components(G, {CC=List()})={
	my(n, s1, s2, s2c, explored, e);
	if(#G[1]==0 || #G[2]==0, return(CC));
	my(data=vector(6));
	rgraph_dfs(~G,~data);
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
	aswell as a ribbon graph with the remaining 
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
			/*Build ribbon graph out of 
			the connected component*/
			listput(~s2CCc, c); 
			foreach(c, j,
				s1CC[j]=s1[j];
			);
		,/*else*/
			/*Build ribbon graph out of 
			the remaining
			connected components*/
			listput(~s2leftc, c);
			foreach(c, j,
				s1left[j]=s1[j];
			);
		);
	);
	my(s2CC, s2left);
	s2CC=cycles_to_perm(n, s2CCc);
	s2left=cycles_to_perm(n, s2leftc);

	listput(~CC, perm_normalize_wrt([s1CC,s2CC],1));

	my(Gleft);
	Gleft=perm_normalize_wrt([s1left,s2left],1);

	/*NOTE: Recursive call on the ribbon graph left*/
	CC=rgraph_connected_components(Gleft, CC);

	return(CC);
}


rgraph_one_face_reduction(G, {withGone=0})={
		my(n,s1,s2);
		[s1,s2]=G;
		n=#s1;
		/*Return if it already has one face.*/
		if(#permcycles(s2)==1, return(G));
		my(data=vector(6), T);
		rgraph_dfs(~G,~data);
		T=data[2];
		
		my(in_T, e);
		in_T=vector(n);
		
		foreach(T, e,
			in_T[e[1]]=1;
			in_T[e[2]]=1;
		); 
		
		/*Build s1_ and s2_*/
		my(s1_, s2_, k, s2inv);
		s1_=vectorsmall(n,i,i);
		s2_=vectorsmall(n,i,i);
		s2inv=s2^-1;
		for(i=1, n, 
			if(in_T[i], next);
			s1_[i]=s1[i];
			k=s2inv[i];
			\\k=s2[i];
			while(in_T[k],
				\\k=s2[s1[k]];
				k=s2inv[s1[k]];
			);
			s2_[i]=k;
		);
		s2_=s2_^-1;
		if(withGone, 
			return([s2_, data, perm_normalize_wrt([s1_,s2_], 1, in_T)]);
		,/*else*/
			return([s2_, data]);
		);
}

makeindex(s2_, seed)={
	my(eindex, e, k);
	eindex=vectorsmall(#s2_);
	e=seed;
	k=1;
	until(e==seed,
			eindex[e]=k;
			e=s2_[e];
			k++;
	);
	return(eindex);
}

/*Assumes G is of genus > 0*/
findab(s2_,s1, seed, eindex)={
	my(a, ainv, ainvindex);
	a=seed; ainv=s1[a]; ainvindex=eindex[ainv];
	my(b, binv, binvindex);
	b=s2_[a]; binv=s1[b]; binvindex=eindex[binv];
	while(binvindex<ainvindex && b!=ainv,
		/*here aindex:=eindex[a]<bindex:=eindex[b]<ainvindex by definition*/
		b=s2_[b];
		binv=s1[b];
		binvindex=eindex[binv];
	);
	return([a,b]);
}

/*
   Performs a single cut and paste on the only cycle of s2_,
 representing the oneface reduction of the graph embedding
 [s1,s2].
*/
cut_and_paste_one(s2_, s1, seedslp, eindex)={
	my(n, n_, a,b,vecalpha, vecbeta, vecgamma, vecdelta);
	n=#s1; n_=permorder(s2_);
	[a,b]=findab(s2_, s1, s2_[seedslp], eindex);
	e=seedslp;

	my(w);
	w=vector(n_, u, e=s2_[e]; e);
	e=a;

	my(card);
	card=(eindex[b]-1)-eindex[a];
	vecalpha=vector(card, u, e=s2_[e]; e);
	e=s2_[e];
	
	card=(eindex[s1[a]]-1)-eindex[b];
	vecbeta=vector(card, u, e=s2_[e]; e);
	e=s2_[e];
	
	card=(eindex[s1[b]]-1)-eindex[s1[a]];
	vecgamma=vector(card, u, e=s2_[e]; e);
	e=s2_[e];
	
	vecdelta=vector(n_-eindex[s1[b]], u, e=s2_[e]; e);
	
	
	my(slpc, slpd);
	slpd=slpconcat([vectoslp(vecalpha,n),[[0,b]], vectoslp(vecbeta,n)],n,,1);
	slpc=slpconcat([slpinvert(vectoslp(vecbeta,n),s1,n), slpinvert(vectoslp(vecgamma,n),s1,n), [[0, a]] ],n,,1);

	\\updated_w=gamma beta c d c^-1d^-1 alpha delta
	\\ a and b are used because they are not in the vecs.
	my(updated_w, updated_s2_, updated_eindex);
	updated_w=concat([vecgamma, vecbeta, [a], [b], [s1[a]], [s1[b]], vecalpha, vecdelta]);
	updated_s2_=cycles_to_perm(n, [updated_w]);
	updated_eindex=makeindex(updated_s2_, updated_w[1]);
		
	return([a,b, vecalpha, vecbeta, vecgamma, vecdelta, slpc, slpd, w, updated_w, updated_s2_, updated_eindex]);
}

/*Deals with the genus 0 case of rgraph_buildpres*/
rgraph_genus0pres(f, s2, s1, slpsgammai,slpgis, pointersgi, {testing=0})={
	my(n);
	n=#s1;
	\\Pointers has size 4*g at this point, that is it points on
	\\2*g topological generators and their inverses 
	\\This filters 2*g topological generators and builds
	\\the appropriate relation.
	my(rels, rel);
	rels=List();
	\\Face relation
	rel=vector(f, u, u);
	listput(~rels, rel);

	\\ Makes slps of gammais into a single slp
   	my(slpgammais, pointersgammai);
	[slpgammais, pointersgammai]=slpconcat(slpsgammai, n,\
		vector(#slpsgammai, u, [#slpsgammai[u]]);
	);
	my(fullslp, pointers);
	\\ Add gis and gammais to main slp
	[fullslp, pointers]=slpconcat(
		\\ SLPS
		[slpgis, slpgammais], n,\
		\\ POINTERS
		[pointersgi, pointersgammai]\
	);

	\\ Add loopfaces to main slp
	my(lenslp_gis_gammai);
	lenslp_gis_gammais=#fullslp;
	fullslp=concat([fullslp,\
		concat(vector(f, u, 
			[[n+pointers[u],n+pointers[f+u]],\
			[-(n+pointers[u]), -1],\
			[n+lenslp_gis_gammais+3*u-2, n+lenslp_gis_gammais+3*u-1]]))\
	]);
	\\ Remove pointers of gis and gammais and add those of
	\\ loopfaces.
	pointers=concat([1],vector(f,u, lenslp_gis_gammais+3*u));  
	
	if(testing,
		return([s1, s2, slpsgammai, slpgammais, slpgis, pointersgi,\
		   		fullslp, pointers, rels]);
	,/*else*/
		return([fullslp, pointers, rels]);
	);
}

/*Build the returned slp of generators and loopfaces.*/
rgraph_buildpres(s2_,s1, seedslp, slpsgammai, slpgis, pointersgi, type, {testing=0})={
	my(n, n_, f);
	n=#s1;
	n_=permorder(s2_);
	if(n_==1, error("Genus 0 is handled in another function : this is a bug"));
	f=(n-n_)/2+1;

	my(e);
	e=seedslp;
		
	my(fullslp, pointers, seen, eindex);
	seen=vectorsmall(n);
	\\Each edge e points to e(1).
	eindex=makeindex(s2_,s2_[seedslp]);

	my(vecalpha, vecbeta, vecgamma, vecdelta,\
		slpalpha, slpbeta, slpgamma, slpdelta,\
	   	   	slpc,slpd);
	my(a,b, w, updated_w, updated_s2_, updated_eindex);
	if(type=="oneword",
		pointers=vectorsmall(n_, i, i);
		w=vectorsmall(n_);
		fullslp=vector(n_, u,\
			   	e=s2_[e];\
				if(!seen[s1[e]], seen[e]=1);\
				w[u]=e;
			   	[0, e]);
	,type=="onehandle",/*else if*/
		\\w_one=gamma beta c d c^-1d^-1 alpha delta
		\\ Recover updated pointers of gis before concatening
		\\ the loopfaces to fullslp and pointers.
		\\ In building pointers we wish to point to the 2*g
		\\ generators and the f loopfaces. As each edge
		\\ appears with its inverse in the one word, we
		\\ point only to the first of the two appearing
		\\ following the w_one.
		[a,b, vecalpha, vecbeta,vecgamma, vecdelta, slpc, slpd,
	   	w, updated_w, updated_s2_, updated_eindex]=\
			cut_and_paste_one(s2_,s1,seedslp,eindex);

		slpalpha=vectoslp(vecalpha, n, 0);
		slpbeta=vectoslp(vecbeta, n, 0);
		slpgamma=vectoslp(vecgamma, n, 0);
		slpdelta=vectoslp(vecdelta, n, 0);

		[fullslp, pointers]=slpconcat(
			\\ SLPS
			[slpgamma, slpbeta, slpc, slpd], n,\
			\\ POINTERS
			[vector(#slpgamma, i,
				if(!seen[s1[slpgamma[i][2]]], seen[slpgamma[i][2]]=1); i),\
			vector(#slpbeta, i,
				if(!seen[s1[slpbeta[i][2]]], seen[slpbeta[i][2]]=1); i),\
			[#slpc], [#slpd]]\
		);
		
		\\ Add slpc^-1 and slpd^-1
		fullslp=concat([fullslp,\
			[[-(n+#slpgamma+#slpbeta+#slpc), -1]],\
		   	[[-(n+#slpgamma+#slpbeta+#slpc+#slpd), -1]]]);
		pointers=concat([pointers,[n+#slpgamma+#slpbeta+#slpc+#slpd+1,n+#slpgamma+#slpbeta+#slpc+#slpd+2]]);
		
		[fullslp, pointers]=slpconcat(
			\\ SLPS
			[fullslp,\
	   		slpalpha, slpdelta], n,\
			\\ POINTERS
			[pointers,
			vector(#slpalpha, i,
				if(!seen[s1[slpalpha[i][2]]], seen[slpalpha[i][2]]=1); i),\
			vector(#slpdelta, i,
				if(!seen[s1[slpdelta[i][2]]], seen[slpdelta[i][2]]=1); i)]\
		);
		s2_=updated_s2_;
		eindex=updated_eindex;
		\\ a and b are now c and d respectively for updated_eindex
		\\ This works as only their index in updated_w is used.
		seen[a]=1;
		seen[b]=1;
	);
	\\Pointers has size 4*g at this point, that is it points on
	\\2*g topological generators and their inverses 
	\\This filters 2*g topological generators and builds
	\\the appropriate relation.
	my(rels, rel);
	rels=List();
	[pointers,rel]=buildrel_and_pointers(pointers, fullslp[1][2], s2_, s1, eindex, seen);
	\\Face relation
	rel=concat(rel, vector(f, u, n_/2+u));
	listput(~rels, rel);

	\\ Makes slps of gammais into a single slp
   	my(slpgammais, pointersgammai);
	[slpgammais, pointersgammai]=slpconcat(slpsgammai, n,\
		vector(#slpsgammai, u, [#slpsgammai[u]]);
	);
	\\ Add gis and gammais to main slp
	[fullslp, pointers]=slpconcat(
		\\ SLPS
		[fullslp,\
		slpgis, slpgammais], n,\
		\\ POINTERS
		[pointers,\
		pointersgi, pointersgammai]\
	);

	\\ Add loopfaces to main slp
	my(lenslp_gens_gis_gammai);
	lenslp_gens_gis_gammais=#fullslp;
	fullslp=concat([fullslp,\
		concat(vector(f, u, 
			[[n+pointers[n_/2+u],n+pointers[n_/2+f+u]],\
			[-(n+pointers[n_/2+u]), -1],\
			[n+lenslp_gens_gis_gammais+3*u-2, n+lenslp_gens_gis_gammais+3*u-1]]))\
	]);
	\\ Remove pointers of gis and gammais and add those of
	\\ loopfaces.
	pointers=concat([pointers[1..n_/2],\
			vector(f,u, lenslp_gens_gis_gammais+3*u)]);  
	
	if(testing,
		if(type=="oneword",
			return([slpsgammai, slpgis, pointersgi,\
		   		fullslp, pointers, rels,\
		   		a, b, eindex,\
		   		s1, s2_, w, seen]);
		,/*else*/
			return([slpsgammai, slpgis, pointersgi,\
			   	fullslp, pointers, rels,\
			   	a, b, eindex,\
			   	s1, s2_, w, updated_w, seen,
			   	vecalpha, vecbeta, vecgamma, vecdelta,\
			   	slpalpha, slpbeta, slpgamma, slpdelta, slpc, slpd]);
		);
	,/*else*/
		return([fullslp, pointers, rels]);
	);
}



/*type="oneword","onehandle","geometric"

Given a ribbon graph G with one face, computes a presentation
of the fundamental group of the dual ribbon graph of the given
type.

The output is a tuple [slp, pointers, rels] with format as describe
in rgraph_buildpres.
 */
rgraph_get_presentation(G, {type="oneword"}, {withdfsfG=0}, {testing=0})={
	if(type=="geometric", /*TODO*/ error("Not implemented yet."); return());
	my(Gdual);
	Gdual=rgraph_dual(G);
	my(s1dual, s2dual, s2dualinv, n);
	[s1dual, s2dual]=Gdual;
	s2dualinv=s2dual^-1;
	n=#s2dual;

	my(s2dual_, dfsfGdual, data, T);
	[s2dual_,data]=rgraph_one_face_reduction(Gdual);
	/*#w=O(g)=n-2(f-1)*/
	[T, dfsfGdual]=data[2..3];
	fGdual=data[6];

	my(s2dualc, f);
	s2dualc=permcycles(s2dual);
	f=#s2dualc;
	
	my(makeslpgammai, maxfacesize = 1);
	/*Encodes a face of Gdual starting at seed and*/
	/*went through in reversed orientation.*/
	makeslpgammai=(u-> 
		my(seed, e, k=1);
		if(u==1,
			seed=1;
		,/*else*/
			seed=T[u-1][2];
		);
		my(slpgammai);
		slpgammai=vector(#s2dualc[fGdual[seed]]);
		e=seed;
		/*On commence à s2dual^-1[seed suivant s2dual^-1*/
		/*et termine juste avant. Puis gi+1 va*/
		/*de s2dual[seed] à j suivant s2dual.*/
		/*le mémo c'est gammai commence en seed */
		/*pas gi, et gi contient pas j.*/
		until(e==seed,
			/*Faut tenir compte des générateurs*/
			if(k==1, 
					slpgammai[k]=[0, s1dual[e]];
		    ,\\else
					slpgammai[k]=[k-1+n, s1dual[e]];
				
			);
			e=s2dualinv[e];
			k++;
		);
		maxfacesize=max(maxfacesize, k-1);
		return(slpgammai);
	);
	my(slpsgammai);
	slpsgammai=vector(f,u, makeslpgammai(u));


   	my(slpsgens, seedlp, n_);
	n_=n-2*(f-1);
	
	my(slpgis, pointersgi, approxtotlength=maxfacesize*n);
	/*straight line program to compute loopfaces paths.*/
	slpgis=vector(approxtotlength);
	/*Used to read the slp : pointersgi[fdfsindex]*/
	/*points to the index of the last edge of the path*/
	/*associated to fdfsindex in slpgis.*/
	pointersgi=vector(f);
	slpgis[1]=[-1,0];
	pointersgi[1]=1;

	my(last, estart, eend, flastindex, flastdfsindex, k);
	k=2;
	for(i=2, f,
		/*
			The path associated to gi is used to join
			the first face to the face 
			of index i in the dfs.
		*/
		if(i==2,
			[estart, eend]=[1, T[1][1]];
			last=1;
		,/*else*/
			flastindex=fGdual[T[i-1][1]];
			flastdfsindex=dfsfGdual[flastindex];
			if(flastdfsindex==1,
				[estart, eend]=[1, T[i-1][1]];
			,/*else*/
				[estart, eend]=[T[flastdfsindex-1][2], T[i-1][1]];
			);
			/*Recover index of last instruction in slp*/
			last=pointersgi[flastdfsindex];
		);

		my(e);
		if(s2dual[estart]==eend,
				/*Add Empty path*/
				/*TODO: transformer en ne rien faire et pas incrémenter k.*/
				slpgis[k]=[0, n+last];
				k++;
		,/*else*/
				e=s2dual[estart];
				until(e==eend,
					/*Unique path contained in face*/
					/*going from s2dual[estart] to */
					/*s2dual^-1[eend] in clockwise */
					/*orientation.*/
					slpgis[k]=[n+last, e];
					e=s2dual[e];
					last=k;
					k++;
				);
		);
		pointersgi[i]=k-1;
	);

	slpgis=slpgis[1..(k-1)];

	my(in_T);
	seed=1;
	in_T=vector(n);
	foreach(T, e,
		in_T[e[1]]=1;
		in_T[e[2]]=1;
	);
	/*As G has one face (in particular s2 has no fixed points)*/
	/*Gdual is reduced, that it, for any e s2dual[e]!=e^-1*/
	/*so that we can compute the genus of Gdual.*/
	seedslp=seed;
	my(ret, genus);
	genus=(#permcycles(Gdual[1]*Gdual[2])-#Gdual[1]/2+f-2)\(-2);
	if(genus,
		/*Fails in genus 0*/
		while(in_T[seedslp],
			seedslp=s2dualinv[s1dual[seedslp]];
		);
		ret=rgraph_buildpres(s2dual_,s1dual, seedslp, slpsgammai, slpgis, pointersgi, type, testing);
	,/*else*/
		ret=rgraph_genus0pres(f, s2dual, s1dual, slpsgammai, slpgis, pointersgi, testing);
	);

	if(withdfsfG, return([ret, dfsfGdual]));
	return(ret);
}


/*permrepr = fonction E-> Sd*/
rgraph_from_permrepr(G, d, permrepr)={
	my(s1, n, s0, s0rev, s1rev);
	s1=G[1];
	n=#s1;
	s0=s1*G[2]^-1;

	/*acts as s0 x id_{1,...,d}*/
	s0rev=vectorsmall(n*d);
	for(i=1, n,
		for(j=1, d,
			s0rev[i+j*n]=s0[i]+j*n;
		);
	);

	/*acts as s1 x (rho o s1) on E x {1,...,d}*/
	s1rev=vectorsmall(n*d);
	for(i=1, n,
		for(j=1, d,
			s1rev[i+j*n]=s1[i]+permrepr[s1[i]][j]*n;
		);
	);
	
	return([s1rev, s1*s0rev^-1]);
}
