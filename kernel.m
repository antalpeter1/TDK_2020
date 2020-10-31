function K = kernel(x1,x2,L)
K = exp(-0.5*(x1-x2)'*L*(x1-x2));