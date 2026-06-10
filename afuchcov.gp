/**/


/*Returns h as slp is not normalized. I.e. it is an slp from*/
/*E(G) to E(Grev) insteadof from elts to E(Grev)*/
afuchcov_raw_spair(X, monodromy)={
	my(G,h, d, Grev, slp, pointers);
	[G,h]=map_from_afuch(X);
	d=#monodromy[1];
	monodromy=vector(#h, i, monodromy[h[i]]);
	[Grev, slp, pointers]=map_raw_spair_from_monodromy(G, d, monodromy);

	return([G,h,d, Grev, slp, pointers]);
}

afuchcov_spair(X, monodromy)={
	my(G,h, d, Gq, slp, pointers, seed);
	[G,h]=map_from_afuch(X);
	d=#monodromy[1];
	m=#monodromy;
	monodromy=vector(#h, i, monodromy[h[i]]);
	[slp, pointers, Gq, seed]=map_spair_from_monodromy(G, d, monodromy);
	\\h is a map from E(G) to elts, slp is from elts to E(Gq)
	slp=slpnormalize(slp, h, m);

	return([G,h,d, Gq, slp, pointers, seed]);
}

\\Given a covering representation of an arithmetic fuchsian group Y,
\\i.e. an afuch X and a monodromy representation associated to the
\\map and side pairing stored in X, returns a one face side pairing obtained directly from the
\\monodromy for Y aswell as a "type" presentation for Y.
afuchcov_presentation(X,monodromy, {type="oneword"},{compose=0}, {eval=0})={
	my(G,h,d,Grev,slpspair, pointersspair,n);
	\\ slpspair maps generators generators associated to side pairing for G
	\\ to generators associated to side pairing for Grev
	[G,h,d, Grev, slpspair, pointersspair]=afuchcov_raw_spair(X,monodromy);
	n=#G[1];
	my(slp, pointers,rels, T,TfGrevdual);
	\\ slp maps generators associated to side pairing for Grev, 
	\\ generating Gamma2, to generators of a presentation of the
	\\ given type.
	[slp, pointers, rels, T, TfGrevdual]=\
			map_topological_presentation(Grev, type);

	if(!compose,
		return([slpspair, pointersspair, slp, pointers, rels]);
	);

	[slp, pointers]=slpcompose([#monodromy, #monodromy], [slp, slpspair], [pointers, pointersspair]);
	my(rels);
	rels=[rels];
	\\slpcompose slpspair et slp
	rels=concat(rels, map_get_ellipticrels(X, map_dual(Grev), h, #pointers, TfGrevdual));
	if(!eval, return([slp, pointers, rels]));

	return([gens, rels]);
}



\\ Given an arithmetic fuchsian group of the form O^x for O 
\\ a maximal order in a quaternion algebra A and a prime pr in
\\ F the center of A. Returns the monodromy representation of
\\ O_0(pr) as a vector of permutations, one for each generator
\\ e in elts where elts=afuchelts(X).
\\WARNING: not deterministic 
afuch_monodromy_from_pr(X, pr)={
	my(A, elts, modpr);
	A=afuchalg(X);
	elts=afuchelts(X);
	modpr=algmodprinit(A,pr);

	\\ [x,y] dans P^1 -> [ax+by, cx+dy]
	\\ a chaque g dans elts -> permutation s(g) dans S_{q+1}
	\\ agir avec g sur chaque elt de P^1 -> [a,1] ou [1, 0]
	\\ g dans elts =(a b) (c d)
	my(F, d, zero, one);
	F=algcenter(A);
	d=poldegree(F.pol)*4;
	zero=vector(d)~;
	one=vector(d, i, if(i==1, 1))~;

	my(ffzero, ffone, ffg, ffmod, v);
	ffzero=algmodpr(A, zero, modpr)[1,1];
	ffone=algmodpr(A, one, modpr)[1,1];
	ffg=ffgen(ffone);
	ffmod=ffone.mod;
	v=variable(ffmod);

	my(g, gmodpr, q,f ffj, gffj, s, monodromy, data=vector(3), P1oo, gj);
	p=modpr[1][1];
	f=poldegree(ffmod); \\residue degree
	data[1]=vector(f); \\represent an element of Fq as tuple of elements of Fp
	data[2]=p;
	data[3]=1;
	q=p^f;
	monodromy=vector(#elts, i, vectorsmall(q+1));
	P1oo=[[ffone, ffzero], ffone, ffzero];
	for(i=1, #elts,
		g=elts[i];

		\\gmodprs=vector(#modprs, u, algmodpr(A, g, modprs[u]));
		gmodpr=algmodpr(A, g, modpr); 
		s=vectorsmall(q+1);	
		for(j=0, q-1,
			\\data[1]==j 
			\\convert j to finite field element 
			ffj=inttoff(data[1], ffg, modpr);
			elP1=[[ffj, ffone],ffone, ffzero];
			\\act with g 
			gelP1=P1homography_act(gmodpr, elP1);
			if(P1isoo(gelP1),
				\\ g.j
				gj=q;
			,/*else*/
				\\ g.j
				gj=fftoint(P1reduce(gelP1)[1][1]);
			);
			s[j+1]=gj+1;
			\\increment j 
			incrbasep(~data);
		);
		gelP1=P1reduce(P1homography_act(gmodpr, P1oo));
		if(P1isoo(gelP1),
			s[q+1]=q+1;
		,/*else*/
			s[q+1]=fftoint(gelP1[1][1])+1;
		);
		\\ FAIRE OO via [ffone, ffzero], implémenter P^1
		monodromy[i]=s;
	);
	return(monodromy);
}

\\ Given an arithmetic fuchsian group associated to
\\ an order O and a level N computes the monodromy
\\ representation of the subgroup O_0(N) 
\\ Currently assumes each prime of the level N in F
\\ splits A and the level
\\ is square free in F.
afuch_congruence(X, N, {flag=0})={
	my(A, F, deg, dA, dF, elts);
	A=afuchalg(X);
	F=algcenter(A);
	deg=poldegree(F.pol);
	dF=F.disc;
	dA=algdisc(A);
	elts=afuchelts(X);

	my(prs);
	if(flag,
		prs=apply((v)->(v[1]), Vec(factor(N)~));
		,/*else*/
		prs=apply((v)->(v[1]), Vec(idealfactor(F, N)~));
	);
	my(modprs);
	modprs=vector(#prs, i, algmodprinit(A, prs[i]));
	monodromy=afuch_congruence_get_monodromy(A, elts, modprs[1]);
	return(monodromy);
}
