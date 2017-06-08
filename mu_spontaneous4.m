function y=mu_spontaneous4(n)
u=rand(n,1);
y=zeros(size(u));

y(u<=0.1)=-1;

idx=find((u>0.9)&(u<=0.95));
if ~isempty(idx)
    y(idx)=mu_factor50(length(idx));
end

idx=find(u>0.95);
if ~isempty(idx)
    y(idx)=-mu_factor50(length(idx));
end