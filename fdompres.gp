/*This file implements various function to retrieve one-word*/
/*and one-handle presentations of any fuchsian group as output*/
/*from the function afuchinit coming from James Rickards's */
/*package Fundamental-domains-for-Shimura-curves available at

 	https://github.com/JamesRickards-Canada/Fundamental-domains-for-Shimura-curves

It also features utility functions to store and retrieve
a Fuchsian group from a file. */





/*Returns a map (or graph embedding) associated*/
/*to the fundamental domain stored in X. In the format */
/*of map.gp*/
map_from_afuch(X,{with_h=1})={
		if(with_h, 
				return(map_from_invol(afuchspair(X)));
		,/*else*/
				return(map_from_invol(afuchspair(X))[1]);
		);
}


modpr_incenter(xmodpr)={
	return(!xmodpr[1,2] && !xmodpr[2,1] && xmodpr[1,1]==xmodpr[2,2]);
}

elliptic_order(A, x, {data=0})={
	my(modprs, maxord);
	if(!data,
		my(p, F, dA, dF, deg, splitprs);
		data=vector(2);
		F=algcenter(A);
		deg=poldegree(F.pol);
		dF=F.disc;
		dA=algdisc(A);

		forprime(tempp=3,,
			if(gcd(dF, tempp)!=1 || gcd(dA, tempp)!=1, next);
			p=tempp;
			break
		);
		splitprs=idealprimedec(F,p);
		modprs = \
			vector(#splitprs, i, algmodprinit(A, splitprs[i]));
		data[1]=modprs;
		/*Worst case bound*/
		data[2]=vecprod(apply((l)->(l/(l-1)), primes(2*deg+1)))*2*deg;
	);
	[modprs, maxord]=data;

	my(ords, xpowmodpr, xmodpr);
	ords=vector(#modprs,i,1);

	for(i=1, #modprs,
		xmodpr=algmodpr(A,x,modprs[i]);
		xpowmodpr=xmodpr;
		while(!modpr_incenter(xpowmodpr) && ords[i]<=maxord,
			xpowmodpr=xpowmodpr*xmodpr;
			ords[i]++;
		);
		if(ords[i]==maxord+1, ords[i]=0);
	);
	ord=lcm(ords);
	return([ord,data]);
}
map_get_ellipticrels(X, Gdual, h, m, dfsfGdual)={
		my(sc, elts, n);
		sc=permcycles(Gdual[2]);
		n=#Gdual[2];

		my(is_ell, ellrels, k=1);
		is_ell=vector(n);
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
		my(elts,A, data);
		elts=afuchelts(X);
		A=afuchalg(X);
		for(i=1, n, 
			if(is_ell[i],
				[is_ell[i],data]=elliptic_order(A, elts[h[i]],data);
				k++;
			);
		);

		my(g);
		g=(m-#sc)/2;
		for(i=1, #sc, 
			c=sc[i];
			if(#c!=1, next);
			listput(~ellrels, [2*g+dfsfGdual[i], is_ell[c[1]]]);
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
		[G,h]=map_from_afuch(X);

		my(ret);
		[ret, dfsfGdual]=map_get_presentation(G,type,1,0);
		/*Add elliptic relations and renormalize the slp*/
		ret[1]=slp_normalize(ret[1],h,#afuchspair(X));
		ret[3]=concat(Vec(ret[3]),map_get_ellipticrels(X,map_dual(G), h, #ret[2],dfsfGdual));
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


afuch_spair_from_monodromy()={
	return();
}
afuch_from_monodromy(X, monodromy)={
	my();
	
	return();
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

afuchstore(X, {storagepath})={
	if(!storagepath,
		storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];
	);

	my(n,sig);
	sig=afuchsignature(X);
	n=fileopen(concat([storagepath,"/", fromsig(sig)]),"a");

	filewrite(n,Str(X));
	fileclose(n);
	return();
}

afuchfromfile(sig,{storagepath})={
	if(!storagepath,
		storagepath=externstr("find /home -wholename */*Cohomology-of-Shimura-curves/storage/fdom")[1];
	);
	n=fileopen(concat([storagepath, "/", tosig(sig)]));
	eval(concat("X= ",fileread(n)));
	return(X);
}

infoafuch(X)={
	my(G);
	G=map_dual(map_from_afuch(X)[1]);
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

