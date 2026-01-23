test_vectoslp({n=25})={
	my(testlin1, testlin2, testpow3, testpow4,\
				v1, v2, v3, v4);

	testlin1=[[0, 3], [n+1, 2], [n+2, 10], [n+3, 5]];
	v1=[3, 2, 10, 5];

	my(slplin);
	slplin1=vectoslp(v1,n);
	if(slplin1!=testlin1, 
		error("test_vectoslp failed.");
	);

	testlin2=[[0, 20], [n+1, 7], [n+2, 15], [n+3, 9]];
	v2=[20, 7, 15, 9];

	my(slplin2);
	slplin2=vectoslp(v2,n);
	if(slplin2!=testlin2, 
		error("test_vectoslp failed.");
	);

	testpow3=[[-25, -30], [n+1, 4], [n+2, 13], [n+3, 7]];
	v3=concat(vector(30, u, -25),[4,13,7]);

	my(slppow3);
	slppow3=vectoslp(v3,n);
	if(slppow3!=testpow3, 
		error("test_vectoslp failed.");
	);

	testpow4=[[-19, 22], [n+1, 21], [n+2, 23], [n+3, 25]];
	v4=concat(vector(22, u, 19),[21, 23, 25]);

	my(slppow4);
	slppow4=vectoslp(v4,n);
	if(slppow4!=testpow4, 
		error("test_vectoslp failed.");
	);
	return();
}

test_slptovec({n=25})={
	my(testlin1, testlin2, testlin3, testlin4,\
			testpow1, testpow2, testpow3, testpow4,\
				v1, v2, v3, v4);

	testlin1=[[0,3], [n+1, 3], [n+2, 2], [n+3, 10], [n+4, 5]];
	testpow1=[[-3, 2], [n+1, 2], [n+2, 10], [n+3, 5]];
	v1=[3, 3, 2, 10, 5];

	my(vlin1,vpow1);
	vlin1=slptovec(testlin1,n);
	vpow1=slptovec(testpow1,n);
	if(vlin1!=v1 || vpow1!=v1, 
		error("test_slptovec failed.");
	);

	testlin2=[[0,20], [n+1, 20], [n+2, 7], [n+3, 15], [n+4, 9]];
	testpow2=[[-20, 2], [n+1, 7], [n+2, 15], [n+3, 9]];
	v2=[20, 20, 7, 15, 9];

	my(vlin2,vpow2);
	vlin2=slptovec(testlin2,n);
	vpow2=slptovec(testpow2,n);
	if(vlin2!=v2 || vpow2!=v2, 
		error("test_slptovec failed.");
	);

	testlin3=concat(vector(30, u, [(u!=1)*n+(u-1), -25]), [[n+30, 4], [n+31, 13], [n+32, 7]]);
	testpow3=[[-25, -30], [n+1, 4], [n+2, 13], [n+3, 7]];
	v3=concat(vector(30, u, -25),[4,13,7]);

	my(vlin3,vpow3);
	vlin3=slptovec(testlin3,n);
	vpow3=slptovec(testpow3,n);
	if(vlin3!=v3 || vpow3!=v3, 
		error("test_slptovec failed.");
	);

	testlin4=concat(vector(22, u, [(u!=1)*n+(u-1), 19]), [[n+22, 21], [n+23, 23], [n+24, 25]]);
	testpow4=[[-19, 22], [n+1, 21], [n+2, 23], [n+3, 25]];
	v4=concat(vector(22, u, 19),[21, 23, 25]);

	my(vlin4,vpow4);
	vlin4=slptovec(testlin4,n);
	vpow4=slptovec(testpow4,n);
	if(vlin4!=v4 || vpow4!=v4, 
		error("test_slptovec failed.");
	);
	return();
}

test_slpconcat({n=25})={
	my(testlin1, testlin2, pointerslin1, pointerslin2,\
		testpow1, testpow2, pointerspow1, pointerspow2,\
		testnonlin1, testnonlin2, pointersnonlin1, pointersnonlin2);
	testlin1=[[0,3], [n+1, 3], [n+2, 2], [n+3, 10], [n+4, 5]];
	pointerslin1=[5];

	testlin2=concat(vector(30, u, [(u!=1)*n+(u-1), -25]), [[n+30, 4], [n+31, 13], [n+32, 7]]);
	pointerslin2=[33];

	testpow1=[[-3, 2], [n+1, 2], [n+2, 10], [n+3, 5],[n+4, 18]];
	pointerspow1=[5];

	testpow2=[[-25, -30], [n+1, 4], [n+2, 13], [n+3, 7]];
	pointerspow2=[4];

	testnonlin1=[[0, 6], [n+1, 1], [-(n+2), 8], [n+3, 10], [n+1, 24], [n+5, 21]];
	pointersnonlin1=[4, 6];

	testnonlin2=[[-23, 4], [-3, 4], [0,20], [n+3, 16], [n+4, 1], [-(n+4), 40]];
	pointersnonlin2=[1, 2 , 5, 6];

	my(slp, pointers);
	[slp, pointers]=slpconcat([testlin1, testlin2], n, [pointerslin1, pointerslin2]);
	realslp=concat([testlin1, vector(30, u, [(u!=1)*(n+5)+(u-1), -25]), [[n+35, 4], [n+36, 13], [n+37, 7]]]);
	realpointers=[5,38];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);

	[slp, pointers]=slpconcat([testlin2, testpow1], n, [pointerslin2, pointerspow1]);
	realslp=concat(testlin2, [[-3, 2], [n+34, 2],[n+35, 10], [n+36, 5], [n+37, 18]]);
	realpointers=[33, 38];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);

	[slp, pointers]=slpconcat([testlin1, testnonlin1], n, [pointerslin1, pointersnonlin1]);
	realslp=concat(testlin1, [[0, 6], [n+5+1, 1], [-(n+5+2), 8], [n+5+3, 10], [n+5+1, 24], [n+5+5, 21]]);
	realpointers=[5, 9, 11];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);
	[slp, pointers]=slpconcat([testpow2, testnonlin2], n, [pointerspow2, pointersnonlin2]);
	realslp=concat(testpow2, [[-23, 4], [-3, 4], [0,20], [n+4+3, 16], [n+4+4, 1], [-(n+4+4), 40]]);
	realpointers=[4, 5, 6, 9, 10];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);

	/*Concat and connect*/
	[slp, pointers]=slpconcat([testpow2, testnonlin2], n, [pointerspow2, pointersnonlin2], 1);
	realslp=concat(testpow2, [[-23, 4], [n+4, n+5], [-3, 4], [0,20], [n+5+3, 16], [n+5+4, 1], [-(n+5+4), 40]]);
	realpointers=[4, 6, 7, 10, 11];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);

	[slp, pointers]=slpconcat([testlin2, testpow1], n, [pointerslin2, pointerspow1], 1);
	realslp=concat(testlin2, [[-3, 2], [n+33, n+34], [n+33+2, 2],[n+33+3, 10], [n+33+4, 5], [n+33+5, 18]]);
	realpointers=[33, 39];
	if(slp!=realslp || pointers!= realpointers,
		error("test_slpconcat failed.");
	);

	return();
}

test_vectoslp();
test_slptovec();
test_slpconcat();
