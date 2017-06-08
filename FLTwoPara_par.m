clear

% Defining constant parameters
N=100; 
u1max=0.7;
u2max=0.3;
beta1max=1;
beta1=0.1;
r2N=3e-3;
N0=100;
KmFold=0.1;
K1sFold=3.3;
K1sFoldMin=3.3;
K2sFold=2;
K2K1=1;
Y1sN=1e3;
Y2sN=1e3;
Y1AN=1e3;
d1=u1max*5e-3;
d2=u2max*5e-3;
T0=17; 

N10=60; 
N20=N0-N10; 
mut_rate=0.01; 
comm_type_num=1; 
comm_rep_num=100; 
max_popul=1e4;
t_bin=0.05;
pcs=1e-15;
K_singular=1e3;
s0=10;

t_binnum=int16(T0/t_bin);

% Initialization
rng('shuffle');
comm_struct=struct('B1_L',zeros(max_popul,1),'B2_L',zeros(max_popul,1),'B1_beta1',zeros(max_popul,1),...
    'B1_KsFold',zeros(max_popul,1),'B1_t',zeros(t_binnum,1),'B2_t',zeros(t_binnum,1),...
    'A',zeros(t_binnum,1),'s',zeros(t_binnum,1),'p',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,'comm_type_num',comm_type_num,'pcs',pcs,'N0',N0);

B1_beta1=zeros(max_popul,1);
B1_KsFold=zeros(max_popul,1);
B1_L=zeros(max_popul,1);
B2_L=zeros(max_popul,1);
death_probability=zeros(max_popul,1);
B1_counter=N10;
B2_counter=N20;
B1_L(1:N10)=1;
B2_L(1:N20)=1;
B1_beta1(1:N10)=beta1;
B1_KsFold(1:N10)=K1sFold;
B1_t=zeros(t_binnum,1);
B2_t=zeros(t_binnum,1);

comm_all(1:comm_type_num*comm_rep_num,1)=struct('B1_L',B1_L,'B2_L',B2_L,'B1_beta1',B1_beta1,...
    'B1_KsFold',B1_KsFold,'B1_t',zeros(t_binnum,1),'B2_t',zeros(t_binnum,1),...
    'A',zeros(t_binnum,1),'s',zeros(t_binnum,1),'p',0,'parentnum',0,'rseed',uint32(0));
rseed=randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
for ri=1:comm_type_num*comm_rep_num
    comm_all(ri).rseed=rseed(ri);
end

% the outer loop that has to be executed serially
for n=1:2;%20*N
    folder_name1=['G' num2str(n)];
    if ~exist(folder_name1,'dir')
        mkdir(folder_name1)
    end  
    % the inner look that can be parallelized
    parfor rep=1:comm_type_num*comm_rep_num
        rnodeseed=comm_all(rep).rseed;
        rng(rnodeseed,'twister');
        comm_rep=comm_struct;
        B1_L=comm_all(rep).B1_L;
        B2_L=comm_all(rep).B2_L;
        B1_beta1=comm_all(rep).B1_beta1;
        B1_KsFold=comm_all(rep).B1_KsFold;
        B1_t=zeros(t_binnum,1);
        B2_t=zeros(t_binnum,1);
        B1_LTemp=zeros(max_popul,1);
        B2_LTemp=zeros(max_popul,1);
        A=zeros(t_binnum,1);
        s=zeros(t_binnum,1);
        p=0;
        
        temp1=find(B1_L>pcs);
        temp2=find(B2_L>pcs);
        B1_counter=length(temp1);
        B2_counter=length(temp2);
        u1maxB1=u1max*B1_L(1:B1_counter);
        u2maxB2=u2max*B2_L(1:B2_counter);
        paras=struct('u1maxB1',u1maxB1,'u2maxB2',u2maxB2,'r2N',r2N,'KmFold',...
            KmFold,'K1sFold',B1_KsFold(1:B1_counter),'K2sFold',K2sFold,...
            'K2K1',K2K1,'Y1AN',Y1AN,'Y1sN',Y1sN,'Y2sN',Y2sN);
        fhandle=@(t,y) as_conc_jacobian(t,y,paras);
        options=odeset('Jacobian',fhandle,'RelTol',1e-5);
        [tx,y]=ode23s(@(t,y) as_conc(t,y,paras),[0 t_bin],[0;s0],options);
        if ~isreal(y)
            error('imaginary value')
        end
        AN=y(:,1);
        sN=y(:,2);
        KsFold1M=B1_KsFold(1:B1_counter)*ones(1,length(tx));
        ANM=ones(B1_counter,1)*AN';
        sNM1=ones(B1_counter,1)*sN';
        u1_coef=(ANM/KmFold)./((ANM/KmFold)+(sNM1./KsFold1M)).*sNM1./(sNM1+KsFold1M)+...
            (sNM1./KsFold1M)./((ANM/KmFold)+(sNM1./KsFold1M)).*ANM./(ANM+KmFold);
        u1=trapz(tx,u1_coef,2)*u1max.*(1-B1_beta1(1:B1_counter));
        u2_coef=sN./(sN+K2K1*K2sFold);
        u2=trapz(tx,u2_coef)*u2max;
        
        B1_LTemp(1:B1_counter)=exp(u1).*B1_L((1:B1_counter));
        B2_LTemp(1:B2_counter)=exp(u2).*B2_L(1:B2_counter);
        p=sum(B1_beta1(1:B1_counter).*(B1_LTemp(1:B1_counter)-B1_L(1:B1_counter))./(1-B1_beta1(1:B1_counter)))+p;
        B1_L(1:B1_counter)=B1_LTemp(1:B1_counter);
        B2_L(1:B2_counter)=B2_LTemp(1:B2_counter);
        A(1)=AN(length(tx));
        s(1)=sN(length(tx));
        
        death_probability=rand(max_popul,1);
        B1_L(death_probability<d1*t_bin)=0;
        B1_beta1(death_probability<d1*t_bin)=0;
        B1_t(1)=sum(B1_L);
        
        death_probability=rand(max_popul,1);
        B2_L(death_probability<d2*t_bin)=0;
        B2_t(1)=sum(B2_L);
        
        div_idx=find(B1_L>=2);
        div_length=length(div_idx);
        if div_length>0
            mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
            B1_beta1(B1_counter+1:B1_counter+div_length)=B1_beta1(div_idx).*mut_multiplier;
            mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
            mut_multiplier=max(mut_multiplier,pcs);
            B1_KsFold(B1_counter+1:B1_counter+div_length)=B1_KsFold(div_idx)./mut_multiplier;
            
            mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
            B1_beta1(div_idx)=B1_beta1(div_idx).*mut_multiplier;             
            mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
            mut_multiplier=max(mut_multiplier,pcs);
            B1_KsFold(div_idx)=B1_KsFold(div_idx)./mut_multiplier; 
            
            B1_L(B1_counter+1:B1_counter+div_length)=B1_L(div_idx)/2;
            B1_L(div_idx)=B1_L(div_idx)/2;
            B1_counter=B1_counter+div_length;
            
            B1_beta1(B1_beta1>beta1max)=beta1max;
            B1_KsFold((B1_KsFold<K1sFoldMin)&(B1_KsFold>pcs))=K1sFoldMin;
        end
        
        div_idx=find(B2_L>=2);
        div_length=length(div_idx);
        if div_length>0            
            B2_L(B2_counter+1:B2_counter+div_length)=B2_L(div_idx)/2;
            B2_L(div_idx)=B2_L(div_idx)/2;
            B2_counter=B2_counter+div_length;
        end
        
        for dt=2:t_binnum
            B1_LTemp=zeros(max_popul,1);
            B2_LTemp=zeros(max_popul,1);
            u1maxB1=u1max*B1_L(1:B1_counter);
            u2maxB2=u2max*B2_L(1:B2_counter);
            
            paras=struct('u1maxB1',u1maxB1,'u2maxB2',u2maxB2,'r2N',r2N,'KmFold',KmFold,...
                'K1sFold',B1_KsFold(1:B1_counter),'K2sFold',K2sFold,...
                'K2K1',K2K1,'Y1AN',Y1AN,'Y1sN',Y1sN,'Y2sN',Y2sN);
            fhandle=@(t,y) as_conc_jacobian(t,y,paras);
            options=odeset('Jacobian',fhandle,'RelTol',1e-5);
            [tx,y]=ode23s(@(t,y) as_conc(t,y,paras),[0 t_bin],[A(dt-1);s(dt-1)],options);
            if ~isreal(y)
                error('imaginary value')
            end
            AN=y(:,1);
            sN=y(:,2);
            KsFold1M=B1_KsFold(1:B1_counter)*ones(1,length(tx));
            ANM=ones(B1_counter,1)*AN';
            sNM1=ones(B1_counter,1)*sN';
            u1_coef=(ANM/KmFold)./((ANM/KmFold)+(sNM1./KsFold1M)).*sNM1./(sNM1+KsFold1M)+...
                (sNM1./KsFold1M)./((ANM/KmFold)+(sNM1./KsFold1M)).*ANM./(ANM+KmFold);
            u1=trapz(tx,u1_coef,2)*u1max.*(1-B1_beta1(1:B1_counter));
            u2_coef=sN./(sN+K2K1*K2sFold);
            u2=trapz(tx,u2_coef).*u2max;
            B1_LTemp(1:B1_counter)=exp(u1).*B1_L((1:B1_counter));
            B2_LTemp(1:B2_counter)=exp(u2).*B2_L(1:B2_counter);
            
            p=sum(B1_beta1(1:B1_counter).*(B1_LTemp(1:B1_counter)-B1_L(1:B1_counter))./(1-B1_beta1(1:B1_counter)))+p;
            B1_L(1:B1_counter)=B1_LTemp(1:B1_counter);
            B2_L(1:B2_counter)=B2_LTemp(1:B2_counter);
            A(dt)=AN(length(tx));
            s(dt)=sN(length(tx));
            
            death_probability=rand(max_popul,1);
            B1_L(death_probability<d1*t_bin)=0;
            B1_beta1(death_probability<d1*t_bin)=0;
            B1_t(dt)=sum(B1_L);
            
            death_probability=rand(max_popul,1);
            B2_L(death_probability<d2*t_bin)=0;
            B2_t(dt)=sum(B2_L);
            
            div_idx=find(B1_L>=2);
            div_length=length(div_idx);
            if div_length>0
                mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
                B1_beta1(B1_counter+1:B1_counter+div_length)=B1_beta1(div_idx).*mut_multiplier;
                mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
                mut_multiplier=max(mut_multiplier,pcs);
                B1_KsFold(B1_counter+1:B1_counter+div_length)=B1_KsFold(div_idx)./mut_multiplier;
                
                mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
                B1_beta1(div_idx)=B1_beta1(div_idx).*mut_multiplier;
                mut_multiplier=mu_spontaneous4(div_length).*double(rand(div_length,1)<=mut_rate)+1;
                mut_multiplier=max(mut_multiplier,pcs);
                B1_KsFold(div_idx)=B1_KsFold(div_idx)./mut_multiplier;
                
                B1_L(B1_counter+1:B1_counter+div_length)=B1_L(div_idx)/2;
                B1_L(div_idx)=B1_L(div_idx)/2;
                B1_counter=B1_counter+div_length;
                
                B1_beta1(B1_beta1>beta1max)=beta1max;
                B1_KsFold((B1_KsFold<K1sFoldMin)&(B1_KsFold>pcs))=K1sFoldMin;
            end
            
            div_idx=find(B2_L>=2);
            div_length=length(div_idx);
            if div_length>0
                B2_L(B2_counter+1:B2_counter+div_length)=B2_L(div_idx)/2;
                B2_L(div_idx)=B2_L(div_idx)/2;
                B2_counter=B2_counter+div_length;
              
            end
            
        end
            temp1=find((B1_L>pcs)&(B1_KsFold<K_singular));
            temp2=find((B2_L>pcs));
            B1_counter=length(temp1);
            comm_rep.B1_L(1:B1_counter)=B1_L(temp1);
            B2_counter=length(temp2);
            comm_rep.B2_L(1:B2_counter)=B2_L(temp2);
            
            comm_rep.B1_beta1(1:B1_counter)=B1_beta1(temp1);
            comm_rep.B1_KsFold(1:B1_counter)=B1_KsFold(temp1);
            
            comm_rep.B1_t=B1_t;
            comm_rep.B2_t=B2_t;
            comm_rep.A=A;
            comm_rep.s=s;
            comm_rep.p=p;
            comm_rep.parentnum=comm_all(rep).parentnum;
            comm_rep.rseed=comm_all(rep).rseed;
            
            comm_all(rep)=comm_rep;
    end
    % save the data 
        folder_name2=['G' num2str(n) '/comm_all'];
        if ~exist(folder_name2,'dir')
            mkdir(folder_name2)
        end
        for i=1:comm_type_num
            comm=comm_all((i-1)*comm_rep_num+1:i*comm_rep_num);
            save([folder_name2 '/comm' num2str(i*comm_rep_num)],'comm');
        end
        % Initialize the next cycle in n
        distrng=rng;
        [~,I]=sort([comm_all.p],'descend');
        comm_all_sorted=comm_all(I);
        rep_counter=0;
        gen_counter=0;
        comm_gen(1:comm_type_num*comm_rep_num,1)=comm_struct;
        for i=1:comm_type_num*comm_rep_num
            if rep_counter>=comm_type_num*comm_rep_num
                break
            end
            dil_factor=floor((comm_all_sorted(i).B1_t(t_binnum)+comm_all_sorted(i).B2_t(t_binnum))/N0);
            if dil_factor==0
                continue
            end
            rep_num_temp=min(dil_factor,comm_rep_num);
            comm_all_idx=min(comm_type_num*comm_rep_num,rep_counter+rep_num_temp);
            comm_all(rep_counter+1:comm_all_idx)=comm_distribute1(comm_all_sorted(i),const_struct,comm_struct,dil_factor,rep_counter,i);
            gen_counter=gen_counter+1;
            comm_gen(gen_counter)=comm_all_sorted(i);
            rep_counter=rep_counter+rep_num_temp;
        end
        comm_gen=comm_gen(1:gen_counter);
        
        rseed=randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
        for ri=1:comm_type_num*comm_rep_num
            comm_all(ri).rseed=rseed(ri);
        end
        
        save([folder_name1 '/comm_gen'],'comm_gen');
        save([folder_name1 '/distrng'],'distrng');
        
end
    
    if ~exist('comm_all','dir')
        mkdir('comm_all')
    end
    for i=1:comm_type_num
        comm=comm_all((i-1)*comm_rep_num+1:i*comm_rep_num);
        save(['comm_all/comm' num2str(i*comm_rep_num)],'comm');
    end
