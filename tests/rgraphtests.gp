/*TODO : Faire des exemples vérifiables à la main.*/
/*Les hardcoder dans les tests ensuite.*/

rgraph_not_connected(n)={
	my(G=rand_rgraph(n));
	while(rgraph_is_connected(G), G=rand_rgraph(n));
	info(G);
	return(G);
}

test_connected_components0(n,{m=1})={
	my(G,CC);
	for(j=1,m,
		G=rgraph_not_connected(n);
		CC=Vec(rgraph_connected_components(G));
		print("\n\n Connected components :\n");
		foreach(CC, G,
			rgraph_info(G,1);
		);
	);
	return();
}

/*NOTE :Deux trois utilitaires*/

is_perm(s)={
	my(n, m=1);
	n=#s;
	rep=vector(n);
	for(i=1,n,
		if(m<s[i],m=s[i]);
		rep[s[i]]++;
	);
	if(m<n, print("Maximum of im(s) : ",m));
	for(j=1, n, if(rep[j]>1, print(j," a deux antécédents\n")));
	return();
}


/*Requires p=3 mod 4 so that s has no fixed points. */
tested_graph(p)={
		if(p%4 != 3, error("p=1 mod 4"); return(););
		my(S, T, s, t, G, Gone); 
		S=vector(p+1);  T=vector(p+1);
		for(i=1, p, T[i]=(i)%(p) +1);
		T[p+1]=p+1;

	   	for(i=1, p-1, S[i]=lift(-Mod(i,p)^-1));
	   	S[p]=p+1;
	   	S[p+1]=p;

		[s,t]=apply(Vecsmall, [S,T]);
		G=[s,t];

		/*Red : p=17 donne un exemple où one_face_reduction renvoie pas une permutation
		  (manque le 4).*/
		if(permorder(s*t)!=3, error("Problème : s*t est pas d'ordre 3."));
		
		my(res);
		res = rgraph_one_face_reduction(G);
		return(res);
}

/*TODO: Rajouter pour test les composantes type ["a","A"] et ["aA"]*/
test_connected_components()={
		my(G, GCC1,GCC2,GCC3);
		G=[Vecsmall([4, 14, 13, 1, 10, 8, 9, 6, 7, 5, 12, 11, 3, 2, 16, 15]), Vecsmall([2, 3, 4, 1, 10, 7, 8, 9, 6, 11, 12, 5, 15, 16, 14, 13])];

		GCC1=[Vecsmall([4, 6, 5, 1, 3, 2, 8, 7]), Vecsmall([2, 3, 4, 1, 7, 8, 6, 5])];
		GCC2=[Vecsmall([2, 1, 4, 3]), Vecsmall([2, 3, 4, 1])];
		GCC3=[Vecsmall([3, 4, 1, 2]), Vecsmall([2, 3, 4, 1])];

		CCsG=rgraph_connected_components(G);
		if(Vec(CCsG)!=[GCC1,GCC2,GCC3], error("rgraph_connected_components failed : invalid connected components."));

		my(Gempty=[Vecsmall([]), Vecsmall([])],CCsGempty);
		CCsGempty=rgraph_connected_components(Gempty); 
		if(CCsGempty!=List([]), error("rgraph_connected_components failed : failed on empty ribbon graph."));
		return();
}
test_covtree(n, {iter=100})={
		my(G,s1,s2,s2c,f,j=0);
		while(j<=iter,
				G=rand_rgraph(n);
				[s1,s2]=G; 
				s2c=permcycles(s2);
				f=#s2c;
				if(f==1,next); 
				if(!rgraph_is_connected(G), next);
				j++;

				my(data);
				data=vector(6);
				rgraph_dfs(G,~data);

				my(covtree,fG); 
				covtree=data[2];
				if(#covtree!=f-1,
						error("Not a covering tree, invalid size.");
						return();
				);
				fG=rgraph_face_index(G);

				my(seen, findex, findexinv);
				seen=vector(f, u, u==data[3][1]);
				for(i=1, #covtree, 
						u=covtree[i][1];
						findex=fG[u];
						findexinv=fG[s1[u]];
						if(findex==findexinv,
							   	error("Invalid edge in tree.");
								return();
						);
						if(seen[findexinv],
								error("Not a tree, found a cycle.");
								return();
						,\\else
								seen[findexinv]=1;
						);
				);
		);
		return();
}
/*TODO: Tester le genre 0 et non réduit.*/
test_one_face_reduction()={
		my(G11tested, G11redtested, covtree11tested,slpdata11tested, G11red,w11, data11);
		G11tested = [Vecsmall([10, 5, 7, 8, 2, 9, 3, 4, 6, 1]), Vecsmall([2, 3, 4, 5, 6, 7, 8, 9, 10, 1])];
		covtree11tested=[[1,10],[2,5],[3,7]];
		G11redtested=[Vecsmall([3, 4, 1, 2]), Vecsmall([2, 3, 4, 1])];

		[w11, data11, G11red]=rgraph_one_face_reduction(rgraph_dual(G11tested),1);

		covtree11=Vec(data11[2]);
		if(covtree11!=covtree11tested, error("rgraph_one_face_reduction failed : Invalid tree."));

		if(G11red!=G11redtested, error("rgraph_one_face_reduction failed : Invalid resulting graph."));
		return();
}

test_get_presentation()={
	my(G, Gdual, data, w, updated_w, a,b, ainv,binv,vecalpha, vecbeta, vecgamma, vecdelta);
	Gdual=[Vecsmall([6, 3, 2, 5, 4, 1, 100, 62, 21, 33, 12, 11, 32, 29, 26, 82, 65, 87, 43, 34, 9, 61, 68, 77, 83, 15, 28, 27, 14, 31, 30, 13, 10, 20, 42, 74, 71, 48, 45, 85, 75, 35, 19, 86, 39, 47, 46, 38, 70, 59, 52, 51, 58, 55, 54, 57, 56, 53, 50, 69, 22, 8, 99, 88, 17, 81, 78, 23, 60, 49, 37, 73, 72, 36, 41, 84, 24, 67, 80, 79, 66, 16, 25, 76, 40, 44, 18, 64, 98, 95, 92, 91, 94, 93, 90, 97, 96, 89, 63, 7, 102, 101]), Vecsmall([101, 6, 3, 2, 5, 4, 1, 100, 62, 21, 33, 12, 11, 32, 29, 26, 82, 65, 87, 43, 34, 9, 61, 68, 77, 83, 15, 28, 27, 14, 31, 30, 13, 10, 20, 42, 74, 71, 48, 45, 85, 75, 35, 19, 86, 39, 47, 46, 38, 70, 59, 52, 51, 58, 55, 54, 57, 56, 53, 50, 69, 22, 8, 99, 88, 17, 81, 78, 23, 60, 49, 37, 73, 72, 36, 41, 84, 24, 67, 80, 79, 66, 16, 25, 76, 40, 44, 18, 64, 98, 95, 92, 91, 94, 93, 90, 97, 96, 89, 63, 7, 102])];
	G=[Vecsmall([6, 3, 2, 5, 4, 1, 100, 62, 21, 33, 12, 11, 32, 29, 26, 82, 65, 87, 43, 34, 9, 61, 68, 77, 83, 15, 28, 27, 14, 31, 30, 13, 10, 20, 42, 74, 71, 48, 45, 85, 75, 35, 19, 86, 39, 47, 46, 38, 70, 59, 52, 51, 58, 55, 54, 57, 56, 53, 50, 69, 22, 8, 99, 88, 17, 81, 78, 23, 60, 49, 37, 73, 72, 36, 41, 84, 24, 67, 80, 79, 66, 16, 25, 76, 40, 44, 18, 64, 98, 95, 92, 91, 94, 93, 90, 97, 96, 89, 63, 7, 102, 101]), Vecsmall([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 1])];
	data=[[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], List([1, 2, 4, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 35, 36, 37, 38, 39, 40, 76, 77, 78, 79, 68, 69, 50, 51, 53, 54, 56, 46, 72, 88, 89, 90, 91, 93, 96, 27, 30, 101]), List([[1, 1], [6, 2, 2], [6, 4, 2], [1, 7], [100, 8, 5], [62, 9, 6], [21, 10, 7], [33, 11, 8], [33, 13, 8], [32, 14, 10], [29, 15, 11], [26, 16, 12], [82, 17, 13], [65, 18, 14], [87, 19, 15], [43, 35, 16], [42, 36, 22], [74, 37, 23], [71, 38, 24], [48, 39, 25], [45, 40, 26], [85, 76, 27], [84, 77, 19], [24, 78, 18], [67, 79, 36], [24, 68, 18], [23, 69, 17], [60, 50, 29], [59, 51, 30], [59, 53, 30], [58, 54, 32], [58, 56, 32], [48, 46, 25], [74, 72, 23], [65, 88, 14], [64, 89, 35], [98, 90, 39], [95, 91, 40], [95, 93, 40], [98, 96, 39], [29, 27, 11], [32, 30, 10], [1, 101]]), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 27, 24, 23, 42, 43, 17, 18, 19, 20, 21, 22, 34, 28, 29, 30, 31, 32, 33, 36, 25, 35, 26, 37, 38, 39, 40, 41, 44], 44, 102];
	w=[63, 22, 34, 83, 66, 99, 44, 75, 49, 86, 25, 61, 70, 81, 41, 20];
updated_w=[44, 75, 49, 86, 25, 34, 83, 66, 63, 22, 99, 61, 70, 81, 41, 20];
	[a,b, ainv, binv]=[63,22,99,61];
	vecalpha=[];
	vecbeta=[22,34,83];
	vecgamma=[44,75,49,86,25];
	vecdelta=[70, 81, 41, 20];

	my(slpc,slpd);

	return();
}



/*TODO: Étendre en test_dfs*/
test_connected_components(); 
test_covtree(25);
test_one_face_reduction();
