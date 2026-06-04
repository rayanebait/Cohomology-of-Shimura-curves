afuchcov_spair_testdata(X, p)={
	my(A, dA, F, pr, f, q, d, monodromy, G, h, Grev, slp, pointers);
	A=afuchalg(X);
	dA=algdisc(A);
	if(gcd(dA, p)!=1, error("Invalid prime in afuchcov_spair_testdata"));

	F=algcenter(A);
	pr=idealprimedec(F, p)[1];
	f=pr[4];
	q=p^f;
	d=q+1;
	monodromy=afuch_monodromy_from_pr(X, pr);
	[G,h,d, Grev, slp, pointers]=afuchcov_spair(X, monodromy);
	return([G,h,d, Grev,slp, pointers]);
}

afuchcov_spair_test(X,{pmax=13})={
	my(A, dA, elts, p);
	A=afuchalg(X);
	dA=algdisc(A);
	elts=afuchelts(X);
	forprime(pp=3, pmax,
		if(gcd(dA, pp)!=1, next);
		p=pp;
		break;
	);
	my(G,h,d,Grev,slp,pointers);
	[G,h,d,Grev,slp,pointers]=afuchcov_spair_testdata(X, p);
	slp=slpnormalize(slp, h, #elts);

	my(gens, rels, s1rev, s2rev, s0rev);
	[s1rev,s2rev]=Grev;
	s0rev=s2rev^-1*s1rev;
	\\ Pas encore s0rev, les relations sont les cycles de s0.
	\\ Et les cycles de s0 là correspondent aux sommets de G.
	\\ Donc aux faces de Gdual !
	rels=permcycles(s1rev*s0rev*s1rev);
	gens=evalslp([A,algmul,algpow,elts], [slp,pointers]);
	evrels=vector(#rels, i, algmulvec(A, gens, rels[i]));

	my(nbzero, nbnonzero);
	foreach(evrels, ev,
		   	if(algincenter(A, ev), nbzero++,nbnonzero++);
	);

	print("Nombre de non nuls : ", nbnonzero," et nombre de nuls : ", nbzero);
	return();
}
/*TODO: Signature [4,[2,2,2,2,2,2,3],0]-> [4,[(2, 6), (3,1)],0]*/
my(X, pol, F, A, pr, J, Or);

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
F=nfinit(pol);
pr=idealprimedec(F,3)[1];
A=alginit(F, [[pr], [0,1]]);
X=afuchinit(A);

[G,h,d,Grev,slp,pointers]=afuchcov_spair_testdata(X, 7);
afuchcov_spair_test(X, 7);
\\pr=idealprimedec(F, 7)[1];

\\my(monodromy, G,h,d,Grev, slp, pointers);
\\monodromy=afuch_monodromy_from_pr(X, pr);
\\[G, h, d, Grev, slp, pointers]=afuchcov_spair(X, monodromy);



