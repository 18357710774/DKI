function D = geodesic_distance_compute(X, Y)

Z = X * Y';
D = acos(Z);
D = real(D);
D = D-diag(diag(D));
