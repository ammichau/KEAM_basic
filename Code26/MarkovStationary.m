function pbar=MarkovStationary(P)
%P=P';
[ns ms]=size(P);
n=ns;
pbar=zeros(n);
while n>1
n1=n-1;
s=sum(P(n,1:n1));
P(1:n1,n)=P(1:n1,n)/s;
n2=n1;
while n2>0
P(1:n1,n2)=P(1:n1,n2)+P(1:n1,n)*P(n,n2);
n2=n2-1;
end
n=n-1;
end
%recursion
pbar=ones(ns,1);
j=2;
while j<=ns
j1=j-1;
pbar(j,1)=sum(pbar(1:j1).*(P(1:j1,j)));
j=j+1;
end
pbar=pbar/(sum(pbar));