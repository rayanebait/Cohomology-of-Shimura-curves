\\ Build a covering map coming from a monodromy
\\ action.
\\ The new edge set is Erev=E x {1,..., d}, we view it as 
\\ {1,...,#E*d} with lexicographic ordering.
\\ The covering vertex permutation s0rev is s0 x id_{1,...,d}
\\ while s1rev is (i,j)->s1(i) x (monodromy(s1(i))(j)).
map_from_monodromy(G, d, monodromy)={
	my(s0, s1, s2, n);
	[s1,s2]=G;
	s0=(s2^-1)*s1;
	n=#s0;

	my(s0rev, s1rev, m=d*n);
	s0rev=vectorsmall(m);
	s1rev=vectorsmall(m);
	for(i=1, d,
		for(j=1, n,
			s0rev[j+(i-1)*n]=s0[j]+(i-1)*n;
			s1rev[j+(i-1)*n]=s1[j]+(monodromy[s1[j]][i]-1)*n;
		);
	);
	return([s1rev, s1rev*s0rev^-1]);
}

\\ Projection map E(Grev)=Ex{1,...,d}->E(G)=E
map_proj(n, d)={
	my(proj);
	proj=vectorsmall(n*d);
	for(i=1, d,
		for(j=1, n,
			proj[j+(i-1)*n]=j;
		);
	);
	return(proj);
}



\\ Given a connected map G->S with S=Gamma1\h with orbifold points at vertices of G
\\ and a covering orbifold Srev -> S with Srev=Gamma2\h and Gamma2 a subgroup
\\ Gamma1 coming from monodromy,
\\ computes a d vertex map Grevone -> Srev with orbifold
\\ points at its vertices as well as a side pairing for Gamma2.
map_raw_spair_from_monodromy(G, d, monodromy)={

	my(Gdual, Gdualrev);
	\\Get single vertex map corresponding to unique polygon P
	Gdual=map_dual(G);
	\\Now has d vertices corresponding to conjugate polygons gP
	\\ for gGamma2 in Gamma2\Gamma1 of size d. 

	Gdualrev=map_from_monodromy(Gdual, d, monodromy);

	my(s0rev, s1rev, s2rev);
	[s1rev,s2rev]=Gdualrev;
	s0rev=s2rev^-1*s1rev;

	my(vG, T, TvG, slpspaths, pointersspath, n,f);
	n=#Gdual[2];
	vG=map_face_index([s1rev,s0rev]);
	\\ Spanning tree T in underlying graph of Gdualrev
	[T,TvG]=map_getT(Gdualrev, 0);
	f=#T+1;

	\\if(f!=d, error(""));


	\\ slp from E(Gdual) to {alphai}={gi}\subset Gamma1.
	\\Yields paths
	\\ alphai=slpspaths[pointersspath[i]] such that if P is a fundamental
	\\ polygon for Gamma1 yielding G with side pairing given by 
	\\ elts as in framework.txt. Then if gi is the element of Gamma1
	\\ associated to alphai using elts then gi.P is a fundamental polygon
	\\ for Gamma1 with side pairing given by gi.elts.gi^{-1} and
	\\ the property that the union U_i giP is a single connected 
	\\ fundamental polygon for Gamma2.


	\\ slp from E(Gdualrev) to paths {alphai}
	[slpspaths, pointersspath]=map_liftalongT(Gdualrev, T, TvG);
	my(proj);
	proj=map_proj(n,d);
	\\ slp from E(Gdual) to paths {proj(alphai)}
	slpspaths=slpnormalize(slpspaths, proj, n);

	my(slptemp, pointerstemp);
	\\ identity slp from E(Gdual) to E(Gdual)
	slptemp=vector(n,i,[0,i]);
	pointerstemp=vector(n,i,i);
	\\ slp from E(Gdual) to E(Gdual) U {alphai}
	[slpspaths, pointersspaths]=slpconcat([slptemp, slpspaths], n, [pointerstemp, pointersspath]);

	my(slprev, pointersrev, e);
	\\ slp from E(Gdual)U{alphai} to E(Gdualrev)
	slprev=vector(3*d*n);
	pointersrev=vector(d*n, u, 3*u);
	for(i=1, d,
		for(j=1, n,
			e=j+(i-1)*n;
			\\ Ici je me trompe c'est pas ça 
			\\ je me balade dans G^* donc pas besoin de vG
			\\ utile pour plus tard
			\\slprev[3*(e-1)+1]=[n+vG[e], j];
			\\slprev[3*(e-1)+2]=[n+vG[s1rev[e]],-1];
			\\slprev[3*(e-1)+3]=[n+f+3*(e-1)+1, n+f+3*(e-1)+2];

			\\alphak^{-1} for path alphak from 1 to vG[e]
			\\here vG[e] is almost always i
			slprev[3*(e-1)+1]=[-(n+TvG[vG[e]]), -1];
			\\alphak^{-1}.e for e as seen in G^*, here written as j
			slprev[3*(e-1)+2]=[n+f+3*(e-1)+1, j];
			\\alphak^{-1}.e.alphak
			slprev[3*(e-1)+3]=[n+f+3*(e-1)+2, n+TvG[vG[e]]];
			pointersrev[e]=3*e;
		);
	);
	my(slp, pointers);
	\\ Side pairing from E(Gdual) to E(Gdualrev)
	[slp,pointers]=slpcompose([n, n+f], [slpspaths, slprev], [pointersspaths, pointersrev]);

	return([[Gdualrev[1], Gdualrev[2]^-1*Gdualrev[1]],slp, pointers]);
}

map_spair_from_monodromy(G, d, monodromy)={
	my(Grev, slp, pointers, Gdualrev);
	\\ Slp from E(G) to E(Grev), side pairing for E(Grev) in terms of 
	\\ E(G)
	[Grev, slp, pointers]=map_raw_spair_from_monodromy(G,d, monodromy);
	Gdualrev=[Grev[1], Grev[1]*Grev[2]^-1];

	my(slpq, pointersq, Gq, seed);
	\\ Gq is the one vertex reduction of Gdualrev
	\\ slp from E(Gdualrev) to E(Gq), side pairing for E(Gq) in terms of
	\\ E(Gdualrev).
	[slpq, pointersq, Gq, seed]=map_quotient_spair([Grev[1], Grev[1]*Grev[2]^-1]);

	\\ slp from E(G) to E(Gq), also side pairing for Gq in terms of G
	\\ seed is returned for relations computations in mappres
	[slp, pointers]=slpcompose([n, d*n], [slp, slpq], [pointers, pointersq]);
	return([slp, pointers, Gq, seed]);
}
