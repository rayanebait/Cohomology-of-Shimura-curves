test_rand_kcyc({iter=10})={
	for(n=10,10^2,
		for(k=1, n,
			for(i=1, iter,
				c=rand_kcyc(n,k);
				if(!iscyc(c, k), error("c not a k-cycle in test_rand_kcycle."));
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
			if(!iscyc(c1), error("c1 is not a cycle in test_as_prodoftwocycs"));
			if(!iscyc(c2), error("c2 is not a cycle in test_as_prodoftwocycs"));
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
			c1=rand_kcyc(n, l);
			c2=rand_kcyc(n, m);

			[c1_,c2_]=cyc_lmtokk([c1,c2]);
			if(c1_*c2_!=c1*c2, error("c1_*c2_!=c1*c2 in test_cyc_lmtokk."));
			if(!iscyc(c1_), error("c1_ is not a cycle in test_cyc_lmtokk."));
			if(!iscyc(c2_), error("c2_ is not a cycle in test_cyc_lmtokk."));

			if(!permorder(c1_)==permorder(c2_), error("c1_ and c2_ not of the same length."));
		);
	);
	return();
}

test_inttoff({n=10^3})={
	my(p, q, f, ffg);
	p=nextprime(n);
	f=floor(log(n));
	q=p^f;
	ffg=ffgen(q);
	v=variable(ffg.pol);

	for(i=1,n,
		for(j=1, f,
			u[j]=random(p);
		);
		if(Pol(u)!=inttoff(u, ffg, v).Pol, error("inttoff failed :", Pol(u), " ", inttoff(u,ffg,v)));
	);
	return();
}
test_ffintconv({n=2^5})={
	my(p, q, f, ffg);
	p=nextprime(n);
	f=floor(log(n));
	q=p^f;
	ffg=ffgen(q);
	v=variable(ffg.pol);

	my(digs);
	for(i=0, q-1,
		digs=digits(i,p);
		if(i!=fftoint(inttoff(digs, ffg, v)), error("ff to int and int to ff conversions failed"));
	);
	return();
}
test_incrbasep({p=nextprime(2^4+random(2^4))})={
	my(q, data=vector(3),f);
	f=floor(log(p))+1;
	data[1]=vector(f); \\represent an element of Fq as tuple of elements of Fp
	data[2]=p;
	data[3]=1;
	q=p^f;
	for(i=0, q-1,
		P=Pol(data[1]);
		if(subst(P, variable(P), p)!=i,
			error("test incrbasep failed.");

		);
		incrbasep(~data);
	);
	return();
}

test_rand_kcyc();
test_as_prodoftwocycs();
test_cyc_lmtokk();
test_ffintconv();
test_incrbasep();
