/*This file implements various function to retrieve one-word*/
/*and one-handle presentations of any fuchsian group as output*/
/*from the function afuchinit coming from James Rickards's */
/*package Fundamental-domains-for-Shimura-curves available at

 	https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves

It also features utility functions to store and retrieve
a Fuchsian group from a file. */





/*Returns a ribbon graph (or graph embedding) associated*/
/*to the fundamental domain stored in X. In the format */
/*of rgraph.gp*/
rgraph_from_afuch(X,{with_h=1})={
		if(with_h, 
				return(rgraph_from_invol(afuchspair(X)));
		,/*else*/
				return(rgraph_from_invol(afuchspair(X))[1]);
		);
}

rgraph_get_ellipticrels(G,ellsorder,m)={
		my(Gdual, sc,s1,s2);
		[s1,s2]=G;
		Gdual=rgraph_dual(G);
		sc=permcycles(Gdual[2]);

		my(is_ell, ellrels, k=1);
		is_ell=vector(#s2);
		ellrels=List();
		/*Mark elliptics which corresponds to cycles*/
		/*of length 1 of s0 with probability 1.*/
		foreach(sc, c, if(#c==1, is_ell[c[1]]=1));

		/*
		   The ordering of the exponents in ellsorder
			corresponds to the ordering of appearance 
			of the elliptics in the only cycle of s2 starting
			at 1.
		*/
		for(i=1, #s2, 
			if(is_ell[i],
				is_ell[i]=ellsorder[k];
				k++;
			);
		);

		my(g);
		g=(m-#sc)/2;
		for(i=1, #sc, 
			c=sc[i];
			if(#c!=1, next);
			listput(~ellrels, [2*g+i, is_ell[c[1]]]);
		);
		return(Vec(ellrels));
}
/*Takes as input a fuchsian group X. The type is */
/*either "oneword" or "onehandle". First computes*/
/*a presentation of the fundamental group associated*/
/*to the dual graph embedding G^* associated to*/
/*the graph embedding G coming from the fundamental*/
/*domain stored in X. */


/*As explained in the paper the output is a straight*/
/*line program with a set of pointers and a set of */
/*relations with the following format : 
	-If the S is a set of generators output by evalslp,
	then S[1..2g] is a set of topological generators,
	if type="onehandle" then S[1], S[2] are canonical loops.
	Further S[2g+1,...,2g+f] are the generators corresponding
	to faces of the dual graph embedding associated to X.

  	- if rels=ret[2], then R=rels[1] is a word of integers
such that the product prod_{i=1}^{4g}S[|R[i]|]^{sign(R[i])} prod_{j=1}^f S[R[i]]
evaluates to 1.
	-The other relations Ri = rels[1+i] are vectors of size 2 such that
	S[Ri[1]]^Ri[2]=1.
 */
/*Letting eval=1 evaluates the slp in the underlying algebra of X*/
/*yielding a *type* presentation for X.*/
afuch_presentation(X, {type="oneword"}, {eval=0})={
		my(G,h);
		/*h est utilisé pour traduire le slp provenant*/
		/*de (G,s2) en slp pour (A,elts).*/
		[G,h]=rgraph_from_afuch(X);

		my(ret);
		ret=rgraph_get_presentation(G,type);
		/*Add elliptic relations and renormalize the slp*/
		ret[1]=slp_normalize(ret[1],h,#afuchspair(X));
		ret[3]=concat(Vec(ret[3]),rgraph_get_ellipticrels(G,afuchsignature(X)[2],#ret[2]));
		if(!eval, return(ret));

		my(slp, pointers, rels);
		[slp, pointers, rels]=ret;

		my(A,elts);
		A=afuchalg(X);
		elts=afuchelts(X);

		my(gens);
		gens=evalslp([A,algmul,algpow,elts], [slp,pointers]);
		/*TODO : rajouter les relations des elliptiques*/
		return([gens,rels]);
}

tosig(v)={
	if(type(v)!="t_STR", return(Str(v)));
	return(v);
}

fromsig(sig)={
	my(v);
	v=sig;
	v[2]=Vec(v[2]);
	return(Str(v));
}

/**/
afuchstore(X)={
	my(n, pathtostorage);
	storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];

	my(sig);
	sig=afuchsignature(X);
	n=fileopen(concat([storagepath,"/", fromsig(sig)]),"a");

	filewrite(n,Str(X));
	fileclose(n);
	return();
}

afuchfromfile(sig)={
	my(n);
	storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];
	n=fileopen(concat([storagepath, "/", tosig(sig)]));
	eval(concat("X= ",fileread(n)));
	return(X);
}

infoafuch(X)={
	my(G);
	G=rgraph_dual(rgraph_from_afuch(X)[1]);
	info(G);
	return();
}

algincenter(A,x)={
	return(alginvol(A,x)==x);
}
algeqmodcenter(A,x,y)={
	my(z);
	z=alginv(A,y);
	z=algmul(A,x,z);
	return(algincenter(A,z));
}

