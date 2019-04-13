function [pdf_poiss] = my_poisspdf(x,lam)

ade = x*log(lam) - lam - gammaln(x+1);
pdf_poiss = exp(vpa(ade));

end

