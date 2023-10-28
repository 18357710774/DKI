function [Nvec, indjCell] = C0SplitDataSimple(D, c0)

Dindi = D>c0;
Vs = setsLambda1(Dindi);       
[Nvec, indjCell] = setsReLambda1(Dindi, Vs);