/*Builds (i j) in Sn */
maketij(n,i,j)={
	my(tij);

	tij=vectorsmall(n,u,u);
	tij[j]=i;
	tij[i]=j;
	return(tij);
}
/*Builds (a b c) in Sn */
make3c(n,a,b,c)={
	my(c3);
	c3=vectorsmall(n,u,u);
	c3[a]=b;
	c3[b]=c;
	c3[c]=a;
	return(c3);
}

/*Random involution with k disjoint transpositions 
in Sn*/
rand_invol(n,k)={
	if(k>n\2, print("invalid transposition number\n"); return(0));

	my(s,m, inv);
	s=vector(n, u, u);
	m=n;
	inv = vectorsmall(n,i,i);

	my(o, j, invj, tij);
	for(i=1, k,
		o = random(m)+1;
		j=s[o];
		s=s[^(o)];
		m=m-1;

		o = random(m)+1;
		invj=s[o];
		s=s[^(o)];
		m=m-1;

		tij = maketij(n,j,invj);
		inv=inv*tij;
	);
	return(inv);
}

/*random k-cycle in Sn*/
rand_kcycle(n,k)={
	if(k>n, print("invalid cycle size\n"); return(0););
	my(s, m, veccyc);

	s=vectorsmall(n,u,u);
	m=n;
	veccyc = vectorsmall(k);

	my(o,j);
	for(i=1,k,
		o = random(m)+1;
		j=s[o];
		s=s[^(o)];
		m=m-1;

		veccyc[i]=j;
		);

	my(cyc);
	cyc = vectorsmall(n, u, u);
	for(i=1, k-1,
		cyc[veccyc[i]]=veccyc[i+1];
	);

	cyc[veccyc[k]]=veccyc[1];
	return(cyc);
}

/*Random permutation in Sn*/
rand_perm(n) =
{
  my(g=vectorsmall(n,i,i),j,tmp);
  for(i=1,n-1,
    j = random([i,n]);
    tmp = g[i];
    g[i] = g[j];
    g[j] = tmp;
  );
  return(g);
}


cycles_to_perm(n,cs)={
	/*Assumes the cycles are disjoint, so that
	we can iterate the single cycle procedure.*/
	my(f, s, card);
	f=#cs;
	s=vectorsmall(n,i,i);
	for(i=1,f,
		cardc=#cs[i];
		for(j=1, cardc-1,
			s[cs[i][j]]=cs[i][j+1];
		);
		s[cs[i][cardc]]=cs[i][1];
	);
	return(s);
}

permconj(s, g)={
	my(n=#s);
	my(sconj=vectorsmall(n,i,i));
	for(i=1, n, sconj[g[i]]=g[s[i]]);
	return(sconj);
}

/*
Given a permutation s of {1,..., n} with fixed points
A. Returns the only permutation g of {1,...,n} such that
g maps A to {1,...,#A} and is order preserving on A and
{1,...,n}-A. This way g*s*g^-1[1...(n-#A)] lies in
the corresponding copy of S_{n-#A} in S_n.
*/
perm_normalizer(s, {V=vector(#s,u,1)},{wrt="fixed points"})={
	my(n, m);
	n=#s;
	m=0;
	for(i=1,n, if(s[i]==i && V[i], m=m+1););

	my(g, k);
	g=vectorsmall(n,i,i);
	k=0;
	for(i=1, n,
		if(s[i]==i && V[i], 
				k=k+1;
				g[i]=n-m+k;
		,\\else
				g[i]=i-k;
		);
	);
	return([g,m]);
}

/*Utility function : if S is a vector of permutations in Sn, 
if S[p] has fixed points I, S[p] lies in an intersection of
stabilizers \cap_i in I Stab(i),
returns a conjugating permutation g such that g S[p].g^{-1}
lies in the intersection \cap_{i=n-#I}^n Stab(i), return
g.S g^-1.
 */
perm_normalize_wrt(S, p, {V=vector(#S[p],u,1)}, {with_g=0}, {wrt="fixed points"})={
	my(s, g, m);
	s=S[p];
	[g,m]=perm_normalizer(s, V, wrt);

	my(gs);
	gs=vector(#S);
	for(i=1, #S,
		gs[i]=permconj(S[i],g)[1..(#s-m)];
	);
	if(with_g, return([gs,g]));
	return(gs);
}

/*i -> i+k mod n*/
perm_iplusk(n,k)={
	return(vectorsmall(n,u, (u-1+k)%n+1));
}



/*If G has less than 26 edges, returns a word representing
 its associated surface.*/
rgraph_to_word(G, {edgetoletter=0})={
	my(s1,s2,n);
	[s1,s2]=G;
	n=#s2;

	if(n>50, return([[],[]]));
	if(n<=1, return(""));

	my(alphabet, alphabetinv);

	alphabet=["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r", "s", "t", "u", "w", "x", "y","z"];
	alphabetinv=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R", "S", "T", "U", "W", "X", "Y","Z"];
	
	my(k);
	if(!edgetoletter,
		edgetoletter=vector(n);
		k=1;
	
		/*Associate a letter or and its inverse to each
		edge in {1,...,n}*/
		for(i=1, n, 
			if(edgetoletter[i]!=0, next);
			edgetoletter[i]=alphabet[k];
			edgetoletter[s1[i]]=alphabetinv[k];
			k=k+1;
		);
	);

	my(s2c, f);
	s2c=permcycles(s2);
	f=#s2c;

	my(words, fG, face);
	words=vector(f,u, vector(#s2c[u]));
	fG=rgraph_face_index(G);
	/*Build words following the cycle structure of
	s2_c which is normalized, cf. rgraph_normalize.*/
	for(i=1,f,
		face=s2c[i];
		for(j=1, #face,
			words[i][j]=edgetoletter[face[j]];
		);
	);

	return([apply(strjoin, words), edgetoletter]);
}

rgraph_info0(G, {w=0},{c=1}, {edgetoletter=0})={
	if(#G[2]<=1, print("\nEmpty ribbon graph\n"); return());
	my(s1, s2, red, conn);
	red=rgraph_is_reduced(G);
	if(c, conn=rgraph_is_connected(G));
	[s1,s2]=G;


	my(v,e,f, g);
	[v,e,f,g]=rgraph_numbers(G);

	my(not, vs, es, fs);
	not = ["not ",""];
	vs = [" vertices, "," vertex, "];
	es = [" edges, "," edge, "];
	fs = [" faces"," face"];

	if(c,
		print("The graph is ", not[conn+1],"connected and ",not[red+1],"reduced.\n");
	,/*else*/
		print("The graph is ", not[red+1],"reduced.\n");
	);
	

	if(c,
		if(red * conn, print("It has ", v , vs[(v==1)+1], e, es[(e==1)+1], f,fs[(f==1)+1]," and genus ", g,".\n"),
			print("It has ", v , vs[(v==1)+1], e, es[(e==1)+1], f,fs[(f==1)+1],".\n"),
		);
	,/*else*/
		if(red, 
			print("It has ", v , vs[(v==1)+1], e, es[(e==1)+1], f,fs[(f==1)+1]," and genus ", g,".\n");
		,/*else*/
			print("It has ", v , vs[(v==1)+1], e, es[(e==1)+1], f,fs[(f==1)+1],".\n");
		);
	);
	if(w, print("Word presentation : ",rgraph_to_word(G,edgetoletter)[1]));
	return();
}

fvec(v)={return((i)->(v[i]))}
find_w(s2_,etol)={
	my(w);
	etol=fvec(etol);
	w=Vec(select((c)->(1<#c), permcycles(s2_))[1]);
	w=concat(apply(etol,w));
	return(w);
}
rgraph_info(G,{edgetoletter=0})={
		my(n=#G[1]); 
		if(n<=50,
				rgraph_info0(G,1,edgetoletter);
		,\\else
				rgraph_info0(G);
		);
		return();
}

rgraph_genword(g)={
	if(g<=0, return([[],[]]));
	my(s1g,s2g);
	s1g=concat(vector(g, i, j=(i-1)*4;\
			   	Vecsmall([3+j, 4+j, 1+j, 2+j])));
	s2g=perm_iplusk(4*g, 1);
	return([s1g,s2g]);
}

\\ Assumes G has one vertex
\\ O(nd)
rand_permrepr(G,d)={
	my(s1,s2,s2c,n, seen);
	[s1,s2]=G;
	s2c=permcycles(G[2]);
	c=s2c[1];
	n=#s1;
	seen=vectorsmall(n);
	my(permrepr);
	foreach(c, e,
		if(seen[e], permrepr[e]=permrepr[s1[e]]^-1);
		permrepr[e]=rand_perm(d);
		
	);
	permrepr[n]=vecprod(permrepr[1..n-1])^-1;
	return(permrepr);
}
\\ faire representation
is_permrepr(G, d, perms)={
	my(permrepr);

	return()
}

rand_afuchpermrepr(sig)={return();}
