function dy=as_conc(t,y,paras)

AN=y(1);
sN=y(2);

u1maxB1=paras.u1maxB1; % a column vector
u2maxB2=paras.u2maxB2; % a column vector
r2N=paras.r2N;
KmFold=paras.KmFold;
K1sFold=paras.K1sFold; % a column vector
K2sFold=paras.K2sFold; % a number
K2K1=paras.K2K1;
Y1AN=paras.Y1AN;
Y1sN=paras.Y1sN;
Y2sN=paras.Y2sN;

u1_coef=(AN/KmFold)./((AN/KmFold)+(sN./K1sFold)).*(sN./K1sFold)./((sN./K1sFold)+1)...
    +(sN./K1sFold)./((AN/KmFold)+(sN./K1sFold)).*(AN/KmFold)./((AN/KmFold)+1);
u2_coef=(sN/K2sFold)/((sN/K2sFold)+K2K1);
dAN=-sum(u1maxB1/Y1AN.*u1_coef)+sum(u2maxB2*r2N*u2_coef);
dsN=-sum(u1maxB1/Y1sN.*u1_coef)-sum(u2maxB2/Y2sN*u2_coef);


dy=[dAN;dsN];