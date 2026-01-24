default("parisize", "16G")

/*File with fields and algebras stored as successive :
	pol 
   [a,b]	

 */
fandapath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fieldsandalgebras/fieldsandalgebras")[1];
storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];


/*Load nsamples fuchsian groups stored file fieldsandalgebras*/
/*of storage directory.*/
afuchget_samples(nsamples)={
		my(m);
		m=fileopen(fandapath);
		
		my(Fs,As,Xs);
		Fs=List(); As=List(); Xs=List();
		
		my(data, dataA,i);
		dataF=dataA=i=1;
		while(dataF && dataA && i<=nsamples,
		    	/*if(data="\n")*/
				dataF = fileread(m); dataA = fileread(m);
				if(!#dataF || !#dataA, return(Xs));
		    	listput(~Fs,nfinit(eval(dataF)));
		    	print("Field ", i, " with polynomial ", Fs[i][1], " and discriminant ", nfdisc(Fs[i][1]),".");
		    	listput(~As, alginit(Fs[i],eval(dataA)));
		    	listput(~Xs, afuchinit(As[i],,1));
		    	i++;
		);
		return(Xs);
}

/*Evaluates th.*/
afuchtest_relation(Xs,{prnt=0},{mul=0}, {type="oneword"})={
	if(mul,
		my(nbfailures,sig, cardxs);
		cardxs=#Xs;
		foreach(Xs, X,
			sig=afuchsignature(X);
			\\if(sig[1]<=2, cardxs--; next);
			if(prnt,
				print("Testing one face relation of fuchsian group with signature :", afuchsignature(X));
			);
			fail=afuchtest_relation(X,prnt,0,type);
			nbfailures+=fail;
		);
		print("Testing ", type, " presentation : ", nbfailures, " failures for ",cardxs," samples.");
		return();
	);
	my(X, A, h, elts);
	X=Xs;
	h=rgraph_from_afuch(X)[2];
	A=afuchalg(X);
	elts=afuchelts(X);

	my(S,R, evrel);
	[S,R]=afuch_presentation(X,type,1);
	evrel=algmulvec(A, S, R[1]);

	my(fail);
	if(!algincenter(A,evrel),
			if(prnt,
				print("Evaluation of relation failed :", evrel);
			);
			fail=1;
	,/*else*/
		if(prnt, print(" Ok."););
		fail=0;
	);
	return(fail);
}

/* A list of fields and algebras are stored in 
STORAGEPATH/fieldsandalgebras and can be used to
compute examples for fuchsian groups with given 
signature. Put store=1 to also store each computed
fundamental domain in a file names after its signature.

The signature doesn't identify congruence arithmetic
fuchsian groups uniquely. The name shall be changed 
in upcoming updates.
 */
afuchsamples({nsamples=10},{comp=0},{store=1})={
	my(Xs);
	Xs=List();
	if(comp,
		Xs=afuchget_samples(nsamples);
		if(store,
		   foreach(Xs, X, afuchstore(X));
		);
	,/*else*/
		sigstrs=externstr(concat(["ls ",storagepath," | grep ',' "]));
		for(i=1, min(nsamples,#sigstrs),
			listput(~Xs,afuchfromfile(eval(sigstrs[i])));
		);
	);
	return(Xs);
}

/*Example of use : 
	Run [ret,toeval]=afuchtestdata(X); eval(toeval);
 	rgraph_info(Gdual);
	datadfs[2]
		->Outputs the word representation of
		the dual graph embedding associated to X
		aswell as the covering tree for the one
		face reduction.
*/
afuchtestdata(X, {type="oneword"},{full=0})={
	my(G,h, Gdual, datadfs,n);
	[G,h]=rgraph_from_afuch(X);
	Gdual=rgraph_dual(G);
	n=#G[1];
	my(v,e,f,g);
	[v,e,f,g]=rgraph_numbers(Gdual);

	datadfs=vector(6);
	rgraph_dfs(Gdual, [1,1,1], ~datadfs);


	my(A, elts);
	A=afuchalg(X);
	elts=afuchelts(X);

	my(ret,toeval);
	ret=[G,h,Gdual,n,datadfs,A,elts];

	if(full,
		if(g,
			if(type=="oneword",
				toeval="[G,h,Gdual,n, datadfs, A,elts, slpsgammai, slpgis, pointersgi, fullslp, pointers, rels, a, b, eindex, s1, s2_, w, seen]=ret";
		,/*else*/
				toeval="[G,h,Gdual,n, datadfs, A,elts, slpsgammai, slpgis, pointersgi, fullslp, pointers, rels, a, b, eindex, s1, s2_, w, updated_w, seen, vecalpha, vecbeta, vecgamma, vecdelta, slpalpha, slpbeta, slpgamma, slpdelta, slpc, slpd]=ret";
			);
		,/*else*/
			toeval="[G,h,Gdual,n,datadfs, A,elts, s1,s2,slpsgammai, slpgammais,slpgis, pointersgi, fullslp, pointers,rels]=ret";
		);
		ret=concat(ret, rgraph_get_presentation(G,type,1));
	,/*else*/
		toeval="[G,h,Gdual,n, datadfs, A,elts, fullslp, pointers, rels]=ret;";
		ret=concat(ret,rgraph_get_presentation(G,type));
	);
	return([ret,toeval]);
}


/*
Test if the words representing the underlying surface 
associated to the one word and one handle presentations
are equal in the fuchsian group associated to X modulo
the center of the algebra afuchalg(X).
*/

afuchtest_comparerelation(X)={
	my(G, h);
	[G,h]=rgraph_from_afuch(X);

	my(n,n_,f);
	n=#G[1];
	f=#permcycles(G[1]*G[2]);
	n_=n-2*(f-1);

	my(A,elts);
	A=afuchalg(X);
	elts=afuchelts(X);

	my(Sone,Rone,Sh,Rh);
	[Sone,Rone]=afuch_presentation(X,"oneword",1);
	[Sh,Rh]=afuch_presentation(X,"onehandle",1);

	if(Sone[n_/2+1..#Sone]!=Sh[n_/2+1..#Sh],
			print("Pas les mêmes elliptiques\n");
			return();
	);

	my(S,R, evrel);
	evrelone=algmulvec(A,Sone,Rone[1][1..n_]); 
	evrelh=algmulvec(A,Sh,Rh[1][1..n_]);

	if(!algeqmodcenter(A, evrelone,evrelh),
		print("Les deux évaluations sont différentes modulo le centre.");
	);

	return();
}

/**/
/*TODO: Signature [4,[2,2,2,2,2,2,3],0]-> [4,[(2, 6), (3,1)],0]*/
my(Xs, pol, F, A, pr, J, Or);
Xs=List();

\\ An example with F=Q(z11)^+, D=(1) and level 32.
pol=y^5+y^4-4*y^3-3*y^2+3*y+1;
F=nfinit(pol);
A=alginit(F,[2, [[],[]],[0, 1/2,1/2,1/2,1/2]]);
pr=idealprimedec(F,2)[1];
J=idealpow(F, pr, 2);
Or=algeichlerbasis(A, pr);
listput(~Xs, afuchinit(A, Or));

/* 
An example with F=Q(sqrt(8)), N_F/Q(D)=9 and level (1). Yields 
a congruence arithmetic Fuchsian group Gamma0^D(1) with signature
(1;3). The code then stores the fundamental domain before
testing the oneword and onehandle presentations.

It then computes various testing utilities.
 */

pol=y^2-8;
F=nfinit(pol);
pr=idealprimedec(F,3)[1];
A=alginit(F, [[pr], [0,1]]);
listput(~Xs,afuchinit(A));


/*Storing and testing*/
foreach(Xs, X, afuchstore(X));
afuchtest_relation(Xs,0,1,"oneword");
afuchtest_relation(Xs,0,1,"onehandle");


/*
   Various signatures corresponding to Fuchsian groups that
   can be computed and stored using afuchsamples(nsamples, 1,1).
*/
my(X,sig);
\\sig=[0,[2,2,2,3],0];
\\sig=[0,[2,3,7],0];
\\sig=[1,[2],0];
\\sig=[1,[2,2],0];
\\sig=[1,[3],0];
\\sig=[1, [2, 2, 2, 2, 2, 2, 3, 3], 0];
\\sig=[2,[],0];
sig=[2,[2],0];
\\sig=[2, [3, 3], 0];
\\sig=[3, [], 0];
\\sig=[4, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3], 0];
\\sig=[67,[],0];

X=afuchfromfile(sig);
\\X=Xs[1];

[ret,toeval]=afuchtestdata(X, "onehandle",1);
eval(toeval);

/*
	The following code computes fundamental domains for 50 algebras
  stored in STORAGEPATH/fieldsandalgebras before testing the
  presentations output by afuch_presentation.
*/
\\my(Xs);
\\Xs=afuchsamples(50);
\\afuchtest_relation(Xs,0,1, "oneword");
\\afuchtest_relation(Xs,0,1, "onehandle");

\\afuchfdom_latex(X,"genus67",0);
