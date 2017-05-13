% y=WrappedGauss(x0,mu,sigma,k)

function y=WrappedGauss(x0,mu,sigma,k)

y=zeros(size(x0));
for i=-k:k
    y=y+normpdf(x0+i,mu,sigma);
end

end

