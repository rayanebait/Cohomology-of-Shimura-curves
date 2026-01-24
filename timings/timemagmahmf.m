/* timemagmahmf.m */

/*
  https://magma.maths.usyd.edu.au/magma/handbook/hilbert_modular_forms
*/

SetLogFile("TimeMagmaHMF.log");

function psi(N)
  fa := Factorization(N);
  ps := 1;
  for i := 1 to #fa do
    ps *:= Norm(fa[i,1])^(fa[i,2]-1) * (Norm(fa[i,1])+1);
  end for;
  return ps;
end function;

_<x> := PolynomialRing(Rationals());
F<t> := NumberField(x^3-3*x-1);
ZF := Integers(F);
A := QuaternionAlgebra<F|-3,4*t>;
OO := MaximalOrder(A);

SetStoreModularForms(F, false);

L := PrimesUpTo(100000,F);

i:=1;
while i le #L do
  N := L[i];
  t := Cputime();
  M := HilbertCuspForms(F,N : QuaternionOrder := OO);
  d := Dimension(M : UseFormula := false);
  t := Cputime(t);
  print Norm(N), psi(N), t;
  if Norm(N) lt 1000 then
    i+:=1;
  else
    i := Round(1.1*i)+1;
  end if;
end while;

