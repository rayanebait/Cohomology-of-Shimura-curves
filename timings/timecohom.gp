/* timecohom.gp */

file = "TimeCohom.log"

idealpsi(nf,N) =
{
  my(fa = idealfactor(nf,N));
  ps = 1;
  for(i=1,#fa[,1],
    my(np = idealnorm(nf,fa[i,1]));
    ps *= np^(fa[1,2]-1) * (np+1)
  );
  ps
};

F = nfinit(t^3-3*t-1);
A = alginit(F,[-3,4*t]);

bnd = 100 000;

L = concat([idealprimedec(F,p,floor(log(bnd)/log(p))) | p <- primes([1,bnd])]);
L = vecsort(L, pr -> idealnorm(F,pr));

ii=1;
{while(ii<=#L,
  N = L[ii];
  ord = algeichlerbasis(A,N);
  t = getabstime();
  X = afuchinit(A,ord,1);
  t1 = getabstime()-t;
  if(afuchsignature(X)[1]>0,
    t = getabstime();
    \\P = afuch_presentation(X,"onehandle");
    P = afuch_presentation(X,"oneword");
    t2 = getabstime()-t;
    print([idealnorm(F,N), idealpsi(F,N), strtime(t1), strtime(t2)]);
    write(file, [idealnorm(F,N), idealpsi(F,N), t1, t2]);
  );
  \\print([idealnorm(F,N), idealpsi(F,N), t1, strtime(t1)]);
  if(idealnorm(F,N) < 1000,
    ii++;
  ,\\else
    ii = round(1.1*ii)+1;
  );
)};

