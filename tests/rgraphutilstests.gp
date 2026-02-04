test_rand_kcycle({iter=10})={
	for(n=10,10^2,
		for(k=1, n,
			for(i=1, iter,
				c=rand_kcycle(n,k);
				if(!is_kcycle(c, k), error("c not a k-cycle in test_rand_kcycle."));
			);
		);
	);
	return();
}
test_as_prodoftwocycs({iter=10^2})={
	my(s,n, c1,c2);
	for(n=10, 10^2,
		for(i=1, iter,
			s=rand_perm(n);
			until(permsign(s)==1,
				s=rand_perm(n);
			);
			[c1,c2]=perm_as_prodoftwocycs(s);
			if(c1*c2!=s, error("c1*c2!=s in test_as_prodoftwocycs"));
			if(!is_cycle(c1), error("c1 is not a cycle in test_as_prodoftwocycs"));
			if(!is_cycle(c2), error("c2 is not a cycle in test_as_prodoftwocycs"));
		);
	);
	return();
}

test_cyc_lmtokk({iter=10^2})={
	my(c1,c2, c1_,c2_, l,m, k);
	for(n=10, 10^2,
		for(i=1, iter,
			l=random(n-1)+1;
			m=random(n-1)+1;
			until(!((l-m)%2),
				l=random(n-1)+1;
				m=random(n-1)+1;
			);
			k=abs((l-m))/2;
			c1=rand_kcycle(n, l);
			c2=rand_kcycle(n, m);

			[c1_,c2_]=cyc_lmtokk([c1,c2]);
			if(c1_*c2_!=c1*c2, error("c1_*c2_!=c1*c2 in test_cyc_lmtokk."));
			if(!is_cycle(c1_), error("c1_ is not a cycle in test_cyc_lmtokk."));
			if(!is_cycle(c2_), error("c2_ is not a cycle in test_cyc_lmtokk."));

			if(!permorder(c1_)==permorder(c2_), error("c1_ and c2_ not of the same length."));
		);
	);
	return();
}

test_rand_kcycle();
test_as_prodoftwocycs();
test_cyc_lmtokk();
