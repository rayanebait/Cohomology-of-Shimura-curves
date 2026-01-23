\r rgraph.gp


\\line = p1, p2, +- avec la convention que slope(p1,p2)
sort_wrt_halfplane(pts, l, pm)={
	my(p1,p2,xsort, ysort,in_p1p2_halfplane);
	[p1,p2]=l;
	\\ plus grand au plus petit
	in_p1p2_halfplane=(u)->(-in_halfplane(p1,p2,u,pm));
	pts=select(in_p1p2halfplane, pts);
	xsort=vecsort(pts, 1);
	ysort=vecsort(pts, 2);

	return([xsort, ysort]);
}
find_lines(pts, l)={
	my(ptstemp, postp_l);
	ptstemp=filter_halfplane(pts,l,-1);
	postp_l=postp()
	if(#ptstemp == 0, return([l]));
	my(candidate_l, postp);
	for(i=1, #ptstemp,
		p=ptstemp[i];
		postp=pos_type(l[1],p);
		if(postp)

		ptstemp=filter_halfplane(pts, [l[1],p],1);
		if(#ptstemp==0,
			l_=[l[1],px];
			break
		);
		py=ysort[i];
		ptstemp=filter_halfplane(pts, [l[1],py],-1);
		if(#ptstemp==0,
			l_=[l[1],py];
			break
		);
		if(i==#pts, return(0));
	);
	return(l_);
}
find_thirdline(pts, l1, l2)={

}

half_span(l, m, pm)={
	my(pts);
	pts=vector(m*m);
	for(i=0, m,
		for(j=0, m-1,
			pts[i+m*j]=[i,j];
		);
	);
	pts=filter_halfplane(pts, l, pm);
	return(pts);
}
\\ faut juste énumerer tout les points entiers
quarter_span(l1, l2, m, pm12)={
	[pm1, pm2]=pm12;
	my(pts, i,j);
	pts=vector(m*m);
	for(i=0, m,
		for(j=0, m-1,
			pts[i+m*j]=[i,j];
		);
	);
	pts=filter_halfplane(pts, l1, pm1);
	pts=filter_halfplane(pts, l2, pm2);
	return(pts);
}
\\ Pas clair ou commencer
enumerate_quarter(l1, l2, facesize)={
	my(m, pts);
	quarter_span(l1, l2, ,[-1,1]);
	
}
\\ memo, si j'ai juste l, c'est facile
\\ si j'ai v1 et v2 deux vecteurs directeurs
\\ calculer produit scalaire, si grand, prendre
\\ bcp de multiples, si petit pas bcp.
\\ 
\\ quand produit scalaire proche de 1 (=c)
\\ on peut tjr supp que entre n*v1 et n*v2
\\ y'a n-1 points à plat. J'pense que y'a
\\ Ah quand 0<c, si ||v1||<||v2||, et c
\\ calcule la proj de v1 sur v2 on a
\\ (1/c)*(n-1) points en nv1 et nv2.
\\
\\Il en faut 2+2 pour un (2+2)-gon, 
\\ il en faut 1+2+2+1 pour un 8-gon,

norm(p)={
	return(sqrt(p[1]^2+p[2]^2));
}

\\ Si premier, faire de taille 3
tex_ngoncoord(pts, l, n)={
	my(p1,p2,m);
	[p1,p2]=l_;
	l_=find_secondline(pts, l);
	m=max(norm(p1),norm(p2));
	my(possible_pts);
	if(l_,
		possible_pts=quarter_span(l,l_, m, [-1,1]);
	,/*else*/
		possible_pts=half_span(l, m, -1);
	)
}

tex_polygonsfromrgraph(G)={
	my(s1,s2,s2c,ngons,n,f);
	n=#s2;
	if(n>25, error("Too big"));
	s2c=permcycles(s2);
	f=#s2c;
	ngonssize=vector(f,u,#s2c[u]);
	if(vecmax(ngonsize)>8, error("Faces are too big"));

	my(words, ordering);
	words=rgraph_to_word(G);
	ordering=words;

	my(data,T);
	data=vector(6);
	rgraph_dfs(G,[1,1,0],~data);
	T=data[2];
	\\ Ordonner les words avec 1... n pour changer en coordonnées.
	\\ Si e est dans T, même arête dans le diagramme, construire
	\\ de proche en proche
	\\ 1) faire un n-gon, 2) recoller au bon endroit

}
