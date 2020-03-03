function [fa, fb] = return_f(a_strat,b_strat,X)
Ca=X(1); Cb=X(2); Ta=X(3); Tb=X(4); N=X(5);
% fa
fa = a_strat(1);
if a_strat(2) ~= 0 % if Nutrient sensor on
    fa = fa + a_strat(2) * (N > a_strat(3));
end
if a_strat(4) ~= 0 % if comp tox sensor on
    fa = fa + a_strat(4) * (Tb > a_strat(5));
end
if a_strat(6) ~= 0 % if self tox on
    fa = fa + a_strat(6) * (Ta > a_strat(7));
end
if a_strat(8) ~= 0 % % if QS on
    fa = fa + a_strat(8) * (Ca > a_strat(9));
end
% fb
fb = b_strat(1);
if b_strat(2) ~= 0 % if Nutrient sensor on
    fb = fb + b_strat(2) * (N > b_strat(3));
end
if b_strat(4) ~= 0 % if comp tox sensor on
    fb = fb + b_strat(4) * (Ta > b_strat(5));
end
if b_strat(6) ~= 0 % if self tox on
    fb = fb + b_strat(6) * (Tb > b_strat(7));
end
if b_strat(8) ~= 0 % % if QS on
    fb = fb + b_strat(8) * (Cb > b_strat(9));
end
end
