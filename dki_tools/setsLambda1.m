function Vs = setsLambda1(Dindi)

V_left = 1:size(Dindi,1);
D_tmp = Dindi;
Vs = cell(1,1);

count = 0;
while ~isempty(V_left)
    count = count + 1;
    V_tmp = setLambdaj1(D_tmp);
    V = V_left(V_tmp);
    Vs{count} = V;

    V_left = setdiff(V_left, V);
    D_tmp = Dindi(V_left, V_left);
end

