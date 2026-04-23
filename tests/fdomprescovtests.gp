/*TODO: Signature [4,[2,2,2,2,2,2,3],0]-> [4,[(2, 6), (3,1)],0]*/
my(X, pol, F, A, pr, J, Or);

\\ An example with F=Q(z11)^+, D=(1) and level 32.
\\pol=y^5+y^4-4*y^3-3*y^2+3*y+1;
\\F=nfinit(pol);
\\A=alginit(F,[2, [[],[]],[0, 1/2,1/2,1/2,1/2]]);
\\pr=idealprimedec(F,2)[1];
\\J=idealpow(F, pr, 2);
\\Or=algeichlerbasis(A, pr);
\\listput(~Xs, afuchinit(A, Or));

/* 
An example with F=Q(sqrt(8)), N_F/Q(D)=9 and level (1). Yields 
a congruence arithmetic Fuchsian group Gamma0^D(1) with signature
(1;3). The code then stores the fundamental domain before
testing the oneword and onehandle presentations.

It then computes various testing utilities.
 */

pol=y^2-8;
F=nfinit(pol);
pr=idealprimedec(F,3)[1];
A=alginit(F, [[pr], [0,1]]);
X=afuchinit(A);

pr=idealprimedec(F, 7)[1];

my(monodromy, G,h,d,Grev, slpspair, pointersspair);
monodromy=afuch_cong_monodromy(X, pr);
[G,h,d, Grev, slpspair, pointersspair]=afuchcov_spair(X, monodromy);

