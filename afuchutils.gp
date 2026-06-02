algincenter(A,x)={
	return(alginvol(A,x)==x);
}
algeqmodcenter(A,x,y)={
	my(z);
	z=alginv(A,y);
	z=algmul(A,x,z);
	return(algincenter(A,z));
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
