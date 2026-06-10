afuchcov_raw_spair_test_data(X, p)={
	my(A, dA); 
	A=afuchalg(X);
	dA=algdisc(A);
	if(gcd(dA, p)!=1, error("Invalid prime in afuchcov_raw_spair_test_data"));

	my(F, pr, f, q, d, monodromy);
	F=algcenter(A);
	pr=idealprimedec(F, p)[1];
	f=pr[4];
	q=p^f;
	d=q+1;
	monodromy=afuch_monodromy_from_pr(X, pr);
	my(G, h, Grev, slp, pointers);
	[G,h,d, Grev, slp, pointers]=afuchcov_raw_spair(X, monodromy);
	return([G,h,d, Grev,slp, pointers]);
}


afuchcov_raw_spair_bench(X,{pmax=100})={
	print("\nBenching afuchcov_raw_spair.\n");
	my(A, dA, p);
	A=afuchalg(X);
	dA=algdisc(A);
	forprime(pp=3, pmax,
		if(gcd(dA, pp)!=1, next);
		p=pp;
		my(G,h,d,Grev,slp,pointers);
		[G,h,d,Grev,slp,pointers]=afuchcov_raw_spair_test_data(X, p);
		
		print("Computation of raw side pairing for covering of degree ", d, " : ok");
	);
	\\print("Nombre de non nuls : ", nbnonzero," et nombre de nuls : ", nbzero);
	return();
}


afuchcov_raw_spair_test(X,{pmax=23})={
	print("\nTesting afuchcov_raw_spair.\n");
	my(A, dA, elts, p);
	A=afuchalg(X);
	dA=algdisc(A);
	elts=afuchelts(X);
	forprime(pp=3, pmax,
		if(gcd(dA, pp)!=1, next);
		p=pp;
		my(G,h,d,Grev,slp,pointers);
		[G,h,d,Grev,slp,pointers]=afuchcov_raw_spair_test_data(X, p);

		my(s2dual, n, s2dualc);
		s2dual=map_dual(G)[2];
		n=#s2dual;
		s2dualc=permcycles(s2dual);

		my(v, rels, c);
		v=#s2dualc;
		rels=vector(d*v);
		for(i=1, d,
			\\ Recall that Erev=Ex{1,...,d} and each cycle of s0 corresponding
			\\ to a vertex of G is a relation for X, or an elliptic a word 
			\\representing an elliptic element for X, if X is our
			\\ fuchsian group.
			\\ Let (e,i) be an edge of Grev above e in G. Then (e,i) maps to
			\\ some gi.g_e.gi^-1 in X if e maps to g_e. Further gi depends only on i.
			\\ In particular, the set of elements in Ex{i} verify the same relations
			\\ as those in E.
			for(j=1, v,
				\\
				c=s2dualc[j];
				rels[j+(i-1)*v]=vectorsmall(#c, k, c[k]+(i-1)*n);
			);
		);
		rels=select((c)->(#c!=1), rels);

		my(gens, evrels);
		gens=evalslp([A,algmul,algpow,elts], [slp,pointers]);
		evrels=vector(#rels, i, algmulvec(A, gens, rels[i]));
		
		my(nbzero, nbnonzero);
		foreach(evrels, ev,
			   	if(algincenter(A, ev), nbzero++,nbnonzero++);
		);
		if(nbnonzero!=0,
				print("Computation of raw side pairing for covering of degree ", d, " : failed");
			,/*else*/
				print("Computation of raw side pairing for covering of degree ", d, " : success");
		);
	);

	\\print("Nombre de non nuls : ", nbnonzero," et nombre de nuls : ", nbzero);
	return();
}

\\TODO: vérifier si la numérotation des arêtes de s2q 
\\ est cohérente avec celle calculée dans map_quotientT
\\ honnêtement je pense que pas du tout.
afuchcov_spair_test_data(X, p)={
	my(A, dA);
	A=afuchalg(X);
	dA=algdisc(A);
	if(gcd(dA, p)!=1, error("Invalid prime in afuchcov_spair_test_data"));
	my(F, pr, f, q, d, monodromy);
	F=algcenter(A);
	pr=idealprimedec(F, p)[1];
	f=pr[4];
	q=p^f;
	d=q+1;
	monodromy=afuch_monodromy_from_pr(X, pr);

	my(G, h, Gq, slp, pointers, seed);
	[G,h,d, Gq, slp, pointers,seed]=afuchcov_spair(X, monodromy);
	my(Gone);
	Gone=perm_normalize_wrt(Gq, 1);
	monodromy=vector(#h, i, monodromy[h[i]]);

	my(s0q, s1q, s2q, index, allrels, rels);
	[s1q, s2q]=Gq;
	s0q=s2q^-1*s1q;
		
	index=map_indexcycle(s0q, seed);
	allrels=select((c)->(c[1]!=0), apply((c)->(veccompose(index, c)), permcycles(s2q)));
	rels=apply((c)->(veccompose(index, c)), permcycles0(s2q));
	
	\\print(#permcycles0(s2q), " ", #permcycles(perm_normalize_wrt([s1q, s2q],1)[2]));

	return([G,h,d, Gq,slp, pointers,seed, rels, allrels, index, monodromy, Gone]);
}

afuchcov_spair_test(X,{pmin=3}, {pmax=53})={
	print("\nTesting afuchcov_spair.\n");
	my(A, F,dA, elts, p);
	A=afuchalg(X);
	F=algcenter(A);
	dA=algdisc(A);
	elts=afuchelts(X);
	forprime(pp=pmin, pmax,
		if(gcd(dA, pp)!=1, next);
		if(#idealprimedec(F,pp)!=poldegree(F.pol), next);
		p=pp;
		my(G,h,d,Gq,slp,pointers, seed, rels, allrels);
		[G,h,d,Gq,slp,pointers,seed, rels, allrels]=afuchcov_spair_test_data(X, p);

		gens=evalslp([A,algmul,algpow,elts], [slp,pointers]);
		evrels=vector(#rels, i, algmulvec(A, gens, rels[i]));
		
		my(nbzero, nbnonzero);
		for(i=1, #evrels,
				ev=evrels[i];
			   	if(algincenter(A, ev), 
					if(#rels[i]>2, print(i, " ", rels[i]););
					nbzero++,
					nbnonzero++);
		);

		print("Computation of side pairing for covering of degree ", d,". Fail : ", nbnonzero,", Success : ", nbzero);
	);

	\\print("Nombre de non nuls : ", nbnonzero," et nombre de nuls : ", nbzero);
	return();
}

afuchcov_spair_bench(X,{pmax=100})={
	print("\nBenching afuchcov_spair.\n");
	my(A, dA, elts);
	A=afuchalg(X);
	dA=algdisc(A);
	elts=afuchelts(X);
	forprime(pp=3, pmax,
		if(gcd(dA, pp)!=1, next);
		my(G,h,d,Gq,slp,pointers, p);
		p=pp;
		[G,h,d,Gq,slp,pointers]=afuchcov_spair_test_data(X, p);

		print("Computation of side pairing for covering of degree ", d, " : ok");
	);

	\\print("Nombre de non nuls : ", nbnonzero," et nombre de nuls : ", nbzero);
	return();
}


/*TODO: Signature [4,[2,2,2,2,2,2,3],0]-> [4,[(2, 6), (3,1)],0]*/
my(X, pol, F, A, pr, J, Or, ab);

\\ An example with F=Q(z11)^+, D=(1) and level 32.
\\pol=y^5+y^4-4*y^3-3*y^2+3*y+1;
\\F=nfinit(pol);
\\A=alginit(F,[2, [[],[]],[0, 1/2,1/2,1/2,1/2]]);
\\pr=idealprimedec(F,2)[1];
\\J=idealpow(F, pr, 2);
\\Or=algeichlerbasis(A, pr);
\\listput(~Xs, afuchinit(A, Or));

/* 
An example with F=Q(sqrt(8)), N_F/Q(D)=9 and level (1). Yields 
a congruence arithmetic Fuchsian group Gamma0^D(1) with signature
(1;3).

Then inputing computing any prime pr in F, computes a side 
pairing for the arithmetic Fuchsian group Gamma0^D(pr).


 */
pol=y^2-8;

\\pol=y^2 - 5 /*[2,[]]*/
\\ab=[-3/2*y - 5/2, -2*y - 9];

\\pol=y^2 - 13 /*[1,[2]]*/
\\ab=[-13, -3*y + 4];

F=nfinit(pol);
\\pr=idealprimedec(F, 7)[1];[1, [3]]


\\A=alginit(F, ab);
A=alginit(F, [[pr], [0,1]]);

X=afuchinit(A);

my(G,h,d,Gq,slp,pointers, rels, allrels, index, monodromy,Gone, elts, gens);
[G,h,d,Gq,slp,pointers, seed, rels, allrels, index, monodromy,Gone]=afuchcov_spair_test_data(X, 17);
elts=afuchelts(X);

\\ Let e be an edge of Gq, that is e is in the support of s0q, assume the support
\\ of s0q is ordered from seed using its only cycle, with ordering index. Then
\\ gens[index[e]] pairs s1q[e] to e.
\\ 
gens=evalslp([A, algmul, algpow, elts], [slp,pointers]);
\\afuchcov_raw_spair_test(X,40);
\\afuchcov_raw_spair_bench(X,17);
\\afuchcov_spair_bench(X,17);
afuchcov_spair_test(X,2, 100);
\\afuchcov_spair_bench(X);

\\my(monodromy, G,h,d,Grev, slp, pointers);
\\monodromy=afuch_monodromy_from_pr(X, pr);
\\[G, h, d, Grev, slp, pointers]=afuchcov_spair(X, monodromy);



