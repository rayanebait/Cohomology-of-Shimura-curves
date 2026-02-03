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


cyctoperm0(n, cs)={
	my(s);
	s=cycles_to_perm(n,[cs[1]]);
	for(i=2, #cs,
		mulpermcyc(~s, ~cs[i]);
	);
	return(s);
}

cycles_to_perm(n,cs)={
	/*Assumes the cycles are disjoint, so that
	we can iterate the single cycle procedure.*/
	my(f, s, cardc);
	f=#cs;
	s=vectorsmall(n,i,i);
	for(i=1, f,
		cardc=#cs[i];
		if(!cardc, next);
		for(j=1, cardc-1,
			s[cs[i][j]]=cs[i][j+1];
		);
		s[cs[i][cardc]]=cs[i][1];
	);
	return(s);
}

/*random k-cycle in Sn*/
rand_kcycle(n,k,{conj=1})={
	if(k>n, print("invalid cycle size\n"); return([]););
	if(k==n, return(cycles_to_perm(n,[rand_perm(n)])));
	my(kcycnorm);
	kcycnorm=cycles_to_perm(n, [rand_perm(k)]);
	if(!conj, return(kcycnorm));
	my(g);
	g=rand_perm(n);
	return(permconj(kcycnorm, g));
}

/*Random involution with k disjoint transpositions 
in Sn*/
rand_invol(n,k)={
	if(k>n\2, print("invalid transposition number\n"); return(0));

	my(invol, g, involcycs);
	g=rand_perm(n);

	involcycs = vector(k, i, Vecsmall([2*i-1, 2*i]));
	invol=cycles_to_perm(n, involcycs);
	invol=permconj(invol, g);

	return(invol);
}




/*g.s.g^-1*/
permconj(s, g)={
	my(n=#s);
	my(sconj=vectorsmall(n,i,i));
	for(i=1, n, sconj[g[i]]=g[s[i]]);
	return(sconj);
}


\\ Conjugates c by s : s c s^-1
conjcycperm(~c, ~s)={
	for(i=1, #c,
		c[i]=s[c[i]];
	);
	return();
}


mulcycperm(~c, ~s)={
	tmp=cycles_to_perm(#s, [c])*s;
	for(i=1, #s, s[i]=tmp[i]);
	return();
}

mulpermcyc(~s, ~c)={
	my(sc, tmp);
	sc=#c;
	if(sc==0, return());
	tmp=s[c[1]];
	for(i=1, sc-1, 
		s[c[i]]=s[c[i+1]];
	);
	s[c[sc]]=tmp;
	return();
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




makeindex(s2, seed)={
	my(eindex, e, k);
	eindex=vectorsmall(#s2);
	e=seed;
	k=1;
	until(e==seed,
			eindex[e]=k;
			e=s2_[e];
			k++;
	);
	return(eindex);
}


is_kcycle(s, k)={
	return(#permcycles(s)==#s-k+1);
}

perm_as_prodoftwocycs_small(s,{k})={
	if(permsign(s)==-1, return([]));
	my(n, g, fixedpts);
	n=#s;
	\\bug pour id
	[g,fixedpts]=perm_normalizer(s);
	s=permconj(s, g);
	if(!k, k=n-fixedpts);
	\\s is now of the form (1 2 ... k1)(k1+1 k1+2...k1+k2)(...)
	c=rand_kcycle(n,k,0);
	n=#s;
	error("");
	while(!is_kcycle(s*c,k),
		c=rand_kcycle(n,k,0);
	);
	s=permconj(s, g);
	c=permconj(c, g);
	return([s*c,c^-1]);
}

\\currently slower than perm_as_profodtwocycs_small 
\\because of mulcycperm which is currently very slow.
perm_as_prodoftwocycs(s)={
	my(si, n);
	n=#s;
	si=s;
	my(gi, gis, cyc3, tmp);
	gis=vector(n-3);
	for(i=1, n-3,
	\\We want gi si gi(1)=2 
	\\gi=(i+1 si(i))
		if(si[i]==i,
			gi=Vecsmall([]);
			gis[i]=gi;
			next;
		,si[i]==i+1,/*else if*/
			gi=Vecsmall([]);
		,/*else*/
			gi=Vecsmall([i+1, si[i]]);
			mulcycperm(~gi,~si);
			mulpermcyc(~si,~gi);
		);
		cyc3=Vecsmall([i+2, i+1, i]);
		\\ gsg(i)
		mulcycperm(~cyc3,~si);
		\\Now si= (i+2 i+1 i). gi-1 si-1 gi-1^-1 
		\\ so that si(i)=i
		\\ At this point si is in Stab(1,2, ...,i)\simeq S_{n-i}
		gis[i]=gi;
	);
	\\ At this point si is an even permutation
	\\	in Stab(1,2, ...,n-3)\simeq S_3. I.e. 
	\\ in <(1 2 3)>. 
	my(ci1, ci2);
	\\ ci1*ci2=si;
	ci1=cycles_to_perm(n, [[n-2, si[si[n-2]], si[n-2]]]);
	ci2=ci1;

	my(i,t);
	for(j=1, n-3,
		i=n-2-j;
		\\gi-1=gis[i]
		gi=gis[i];

		\\ci-1,1= gi-1.((i+2 i)((i+1 i).ci1).(i+2 i))gi-1
		t=Vecsmall([i+1,i]);
		mulcycperm(~t, ~ci1);
		t=Vecsmall([i+2, i]);
		mulcycperm(~t, ~ci1);
		mulpermcyc(~ci1, ~t);
		mulpermcyc(~ci1, ~gi);
		mulcycperm(~gi, ~ci1);

		\\ci-1,2=gi-1.((i+2 i)ci,2)gi-1
		t=Vecsmall([i+2, i]);
		mulcycperm(~t, ~ci2);
		mulpermcyc(~ci2, ~gi);
		mulcycperm(~gi, ~ci2);
	);
	
	return([ci1,ci2]);
}

permcycles0(s)={
	return(select((c)->(1<#c), permcycles(s)));
}
cyc_lmtokk(cs)={
	my(ords);
	ords=apply(permorder, cs);

	my(m1, m2);
	m1=ords[1];
	m2=ords[2];
	if(m1==m2, return(cs));
	if(!(m1-m2)%2, return([]));

	my(smallest, m=max(m1,m2));
	if(m==m1, smallest=2, smallest=1);

	my(smallc, bigc, diff=abs(m2-m1), S, n);
	smallc=cs[smallest]; bigc=cs[-smallest%3];
	n=#bigc;
	S=vectorsmall(n);

	\\compute intersection of supports of bigc and smallc
	for(i=1, n,
		if(bigc[i]!=i,
			if(smallc[i]!=i,
				S[i]=1;
			);
		);
	);

	my(cardS, ind);
	cardS=vecsum(S);
	if(cardS==0, 
		\\ The cycles are disjoint
		if(smallest==1,
			smalle=permcycles0(smallc)[1][1];
			bige=permcycles0(bigc)[1][1];
			t=maketij(n, smalle, bige);
		);
	);
	ind=1;
	\\ if S is the intersection of supports of bigc and smallc
	\\ and if i is in S but bigc(i) is not in S, then put
	\\ smallc=smallc.(i bigc(i))
	\\ bigc=(i bigc(i)).bigc
	\\ then i is not in S anymore but bigc(i) is in S
	for(i=1, n,
		if(S[i] && !S[bigc[i]],
			lts[ind]=[i, bigc[i]];
		);
	);

	\\TODO: trouver un bon critère pour le faire sans se soucier.
	my(t,tinv);
	t=vectorsmall(n,i,i);
	for(i=1, #ts,
		mulpermcyc(~t,~ts[i]);
	);
	tinv=t^-1;
	print(permcycles0(t));
	if(smallest=2,
		return([bigc*tinv,tinv*smallc]);
	,/*else*/
		return([smallc*t, t*bigc]);
	);
}

\\TODO : sort 1... n à chaque étape en mettant les sii à la fin,
\\ pour choisir un j pas dans le support, remplacer i+2.
perm_as_prodoftwocycs1(s)={
	my(n, sc, f, g, fixedpts, findex);
	n=#s;
	sc=permcycles(s);
	f=#sc;
	g=vectorsmall(n,i,i);
	findex=vectorsmall(n);
	fixedpts=0;

	my(k);
	k=0;
	for(i=1, #sc,
		c=sc[i];
		if(#c==1, fixedpts++);
		for(j=1, #c,
			g[c[j]]=j+k;
			findex[j+k]=i;
		);
		k+=#c;
	);
	\\s is now of the form (1 2 ... k1)(k1+1 k1+2...k1+k2)(...)
	s=permconj(s, g);
	sc=permcycles(s);

	\\If n-f_n-2 is even, then #cs=1 or 3 otherwise 
	\\ 2 and we should go one more step.
	my(fn_2, fn_1, fn);
	fn_2=f-2*((findex[n-2]==findex[n-1]) && (findex[n-1]!=findex[n])) \
	   	-2*((findex[n-2]!=findex[n-1]) && (findex[n-1]==findex[n])) \
	   	-3*((findex[n-2]!=findex[n-1]) && (findex[n-1]!=findex[n])) \
	   	-((findex[n-2]==findex[n-1]) && (findex[n-1]==findex[n]));
	fn_1=fn_2+(findex[n-2]!=findex[n-1]);
	fn=fn_1+(findex[n-1]!=findex[n]);
		
	\\si=(i si-1(i))si-1 when i isnt a fixed point of si-1, si=si-1 otherwise
	my(sigsn_3,sigsn_2,sigsn_1);
	sigsn_3=(-1)^((n-(3+fn_2))%2)*(permsign(s));
	sigsn_2=(-1)^((n-(2+fn_1))%2)*(permsign(s));
	sigsn_1=(-1)^((n-(1+fn))%2)*(permsign(s));

	\\print("sigsn_3 : ", sigsn_3);
	\\print("sigsn_2 : ", sigsn_2);
	\\print("sigsn_1 : ", sigsn_1);
	\\print("fn_2 : ", fn_2);
	\\print("fn_1 : ", fn_1);
	\\print("fn : ", fn);

	my(c, cs, j, b);
	cs=List();
	c=List();
	b=2*(sigsn_3==1) \
	  +1*(sigsn_3==-1 && sigsn_2==1);
	  
	for(i=0, b, 
		listput(~c, n-b+i);
		j=s[n-b+i];
		if(j!=n-b+(i+1),
			listput(~cs, Vecsmall(c));
			c=List();
		);
	);

	my(i,m, j, lis, ris, f0, f1);
	m=0; j=1; f0=0; f1=f0;
	lis=vector(n-(1+b+fn_2*(sigsn_3==1)\
					+fn_1*(sigsn_3==-1 && sigsn_2==1)\
					+fn*(sigsn_3==-1 && sigsn_2==-1 && sigsn_1==1)));
	ris=vector((fn_2*(sigsn_3==1)\
					+fn_1*(sigsn_3==-1 && sigsn_2==1)\
					+fn*(sigsn_3==-1 && sigsn_2==-1 && sigsn_1==1)));
	\\ m goes through the last index of each cycle
	\\ j goes through every index that is not the last index
	\\ of a cycle -> n-fn-3 in total
	\\ f0 is the index of the cycle we are in in sc
	\\ m+#sc[f0+1]+1 is the first index in the f0+1'th cycle if
	\\ c is the f0'th cycle
	for(i=1, n-b-1,
		f1=findex[i];
		if(f0!=f1,
			f0=f1;
			c=sc[f0];
		);
		if(i==c[#c] && f0<=fn_1,
			m+=#c;
			ris[f0]=Vecsmall([m, sc[f0+1][#sc[f0+1]]]);
		,i!=c[#c],/*else if*/
			lis[j]=Vecsmall([i, i+1]);
			j++;
		);
	);

	\\print(cs);
	c=cycles_to_perm(n, cs);
	my(c01,c02);
	if(sigsn_3==-1,
		c01=maketij(n,n-1,n);
		c02=c01;
	,/*else*/
		if(c[n-2]==n-2,
			c01=cycles_to_perm(n, [[n-2,n-1,n]]);
			c02=c01^-1;
		,/*else*/
			c01=c^2;
			c02=c01;
		);
	);

	\\print("s : ", permcycles(s),"\n");
	\\print("lis : ", lis,"\n");
	\\print("ris : ", ris,"\n");

	my(c1);
	c1=vectorsmall(n,i,i);
	for(j=1, #lis,
		mulpermcyc(~c1, ~lis[j]);
	);
	c1=c1*c01;

	my(c2,j);
	c2=vectorsmall(n,i,i);
	for(l=1, #ris,
		j=#ris-l+1;
		mulpermcyc(~c1, ~ris[l]);
		mulpermcyc(~c2, ~ris[j]);
	);
	c2=c2*c02;

	g=g^-1;
	\\return(cyc_lmtokk(n,c1,c2));
	return(apply((c)->(permconj(c,g)), [c1,c2]));
}



perm_as_commutator(s)={
	my(c1,c2);
	[c1,c2]=perm_as_prodoftwocycs(s);

	return();
}

\\ Assumes G is a one handle word
\\ O(nd)
rand_rgraphpermrepr(c,n,d)={
	my(n, seed, eindex);
	seed=c[1];
	eindex=makeeindex(cycles_to_perm(n, [c]), seed);

	my(eend, eendinv, seen);
	eend=c[#c];
	eendinv=s1[eeinv];
	seen=vectorsmall(n);

	my(permrepr, e);
	permrepr=vector(#c);
	for(i=1, eindex[eendinv],
		e=c[i];
		if(seen[e], next);
		permrepr[e]=rand_perm(d);
		permrepr[s1[e]]=permrepr[e]^-1;
		seen[e]=1;
		seen[s1[e]]=1;
	);
	w1=vecprod(permrepr[seed..eindex[eendinv]-1]);
	w2=permconj(w1, permrepr[eend]);

	for(i=eindex[eendinv]+1, #c-2,
		e=c[i];
		if(seen[e], next);
		permrepr[e]=rand_perm(d);
		permrepr[s1[e]]=permrepr[e]^-1;
		seen[e]=1;
		seen[s1[e]]=1;
	);



	permrepr[n]=vecprod(permrepr[1..n-1])^-1;
	return(permrepr);
}

\\ faire representation
is_permrepr(G, d, perms)={
	my(permrepr);

	return();
}

rand_afuchpermrepr_fromR(R,s1)={

	return();
}
