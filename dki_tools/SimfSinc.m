function y = SimfSinc(x, c1, c2)

r = c1*sqrt(sum(x.^2,2));
y = sin(c2*pi*r)./(r);

end