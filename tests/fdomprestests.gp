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
		print("Testing ", #Xs," samples.");
		foreach(Xs, X,
			sig=afuchsignature(X);
			\\if(sig[1]<=2, cardxs--; next);
			if(prnt,
				print("Testing one face relation of fuchsian group with signature :", afuchsignature(X));
			);
			fail=afuchtest_relation(X,prnt,0,type);
			if(fail, print(sig));
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

	my(S,R, evrel, evells);
	[S,R]=afuch_presentation(X,type,1);
	evrel=algmulvec(A, S, R[1]);
	evells=vector(#R-1, i, algpow(A, S[R[i+1][1]], R[i+1][2]));

	my(fail,ev);
	if(!algincenter(A,evrel),
			if(prnt,
				print("Evaluation of relation failed :", evrel);
			);
			fail=1;
	,/*else*/
		if(prnt, print(" Ok."););
		fail=0;
	);
	for(i=1, #evells,
		/*
		if(algincenter(A, S[R[i][1]]),
			if(prnt,
				print("Relation pointing to trivial generator :", S[R[1]]);
			);
			fail=1;
		,\\else
			fail=0;
		);*/
		ev=evells[i];
		if(!algincenter(A,ev),
			\\print(i, ev); error("");
			if(prnt,
				print("Evaluation of elliptic relations failed :", ev);
			);
			fail=1;
		,/*else*/
			fail=0;
		);
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
	my(Xs,storagepath);
	Xs=List();
	storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];
	if(comp,
		print("Computing ",nsamples," samples.");
		Xs=afuchget_samples(nsamples);
		print("Storing ",nsamples," samples.");
		if(store,
		   foreach(Xs, X, afuchstore(X,storagepath));
		);
	,/*else*/
		sigstrs=externstr(concat(["ls ",storagepath," | grep ',' "]));
		print("Loading ", min(nsamples,#sigstrs)," samples.");
		for(i=1, min(nsamples,#sigstrs),
			listput(~Xs,afuchfromfile(eval(sigstrs[i]),storagepath));
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
afuchtestdata(X, {type="oneword"})={
	my(G,h, Gdual, datadfs,n);
	[G,h]=rgraph_from_afuch(X);
	Gdual=rgraph_dual(G);
	n=#G[1];
	my(v,e,f,g);
	[v,e,f,g]=rgraph_numbers(Gdual);

	my(s2_, data, etol, w, datadfs, Gone);
	[s2_,data,Gone]=rgraph_one_face_reduction(Gdual,1);
	etol=rgraph_to_word(G)[2];
	if(g && #etol,
		w=find_w(s2_,etol);
		wis=rgraph_to_word(Gdual,etol)[1];
	,/*else*/
		w=s2_;
		wis=permcycles(Gdual[2]);
	);

	my(A, elts);
	A=afuchalg(X);
	elts=afuchelts(X);

	my(ret,toeval);
	ret=[G,h,Gdual,n,datadfs,A,elts];

	toeval="[G,h,Gdual,n, datadfs, A,elts, fullslp, pointers, rels,wis,w,etol]=ret;";
	ret=concat(ret, afuch_presentation(X,type));
	ret=concat(ret, [wis,w,etol]);

	
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
(1;3). The code then stores the fundamental domain before
testing the oneword and onehandle presentations.

It then computes various testing utilities.
 */

\\pol=y^2-8;
\\F=nfinit(pol);
\\pr=idealprimedec(F,3)[1];
\\A=alginit(F, [[pr], [0,1]]);
\\listput(~Xs,afuchinit(A));


/*Storing and testing*/
\\foreach(Xs, X, afuchstore(X));
\\afuchtest_relation(Xs,0,1,"oneword");
\\afuchtest_relation(Xs,0,1,"onehandle");
\\
\\
\\/*
\\   Various signatures corresponding to Fuchsian groups that
\\   can be computed and stored using afuchsamples(nsamples, 1,1).
\\*/
my(X,sig);
sig=[0,[2,2,2,3],0];
\\sig=[0,[2,3,7],0];
\\sig=[1,[2],0];
\\sig=[1,[2,2],0];
\\sig=[1,[3],0];
\\\\sig=[1, [2, 2, 2, 2, 2, 2, 3, 3], 0];
\\\\sig=[2,[],0];
\\\\sig=[2,[2],0];
\\sig=[2, [3, 3], 0];
\\sig=[3, [], 0];
\\sig=[4, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3], 0];
\\sig=[12, [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3], 0]
\\sig=[67,[],0];
\\
X=afuchfromfile(sig);
Xs=[X];
\\
\\/*
\\    The following code computes fundamental domains for 
\\    50 algebras stored in STORAGEPATH/fieldsandalgebras before
\\   	testing the presentations output by afuch_presentation.
\\
\\    Running Xs=afuchsamples(50,1,1) first computes the 
\\    fundamental domains before storing them. 
\\    Running Xs=afuchsamples(50,0,1) uses samples stored in
\\    STORAGEPATH/fdom. 
\\    Running Xs=afuchsamples(50,1,0) computes the fundamental
\\    domains without storing them.
\\*/
\\
\\/*Computes 50 fundamental domains*/
\\\\Xs=afuchsamples(50,1,1);
\\
\\/*Retrieves at most 50 fundamental domains from storage*/
\\\\Xs=afuchsamples(50,0,1);
afuchtest_relation(Xs,0,1, "oneword");
\\\\afuchtest_relation(Xs,0,1, "onehandle");
\\
[ret,toeval]=afuchtestdata(X, "oneword");
eval(toeval);
\\\\afuchfdom_latex(X,"genus67",0);
