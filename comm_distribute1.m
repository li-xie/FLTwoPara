function comm_rep=comm_distribute1(comm_all_struct,const_struct,comm_struct,dil_factor,rep_counter,parentnum)
comm_rep_num=const_struct.comm_rep_num;
comm_type_num=const_struct.comm_type_num;
pcs=const_struct.pcs;
% t_binnum=const_struct.t_binnum;
% max_popul=const_struct.max_popul;
% comm_struct=struct('B1_L',zeros(max_popul,1),'B2_L',zeros(max_popul,1),'B1_beta1',zeros(max_popul,1),...
%     'B1_t',zeros(t_binnum,1),'B2_t',zeros(t_binnum,1),'c',zeros(t_binnum,1),'A',zeros(t_binnum,1),'p',0,'parentnum',0,'rseed',uint32(0));
comm_rep(1:comm_rep_num,1)=comm_struct;
B1_L=comm_all_struct.B1_L;
B2_L=comm_all_struct.B2_L;
B1_beta1=comm_all_struct.B1_beta1;
% B1_umax=comm_all_struct.B1_umax;
% B1_KmFold=comm_all_struct.B1_KmFold;
B1_KsFold=comm_all_struct.B1_KsFold;
% B2_umax=comm_all_struct.B2_umax;
% B2_KsFold=comm_all_struct.B2_KsFold;
B1_counter=nnz(comm_all_struct.B1_L>pcs);
B2_counter=nnz(comm_all_struct.B2_L>pcs);
% B1_0=comm_all_struct.B1_t(1);
% B2_0=comm_all_struct.B2_t(1);
% B1_T=comm_all_struct.B1_t(t_binnum);
% B2_T=comm_all_struct.B2_t(t_binnum);

% dil_factor=round(B1_T+B2_T)/(B1_0+B2_0);
rand_temp=ceil(rand(B1_counter+B2_counter,1)*dil_factor);
if dil_factor>=comm_rep_num
    for i=1:comm_rep_num
        temp_idx=find(rand_temp==i);
        B1_num=nnz(temp_idx<=B1_counter);
        B2_num=nnz(temp_idx>B1_counter);
        if B1_num>=1
            comm_rep(i).B1_L(1:B1_num)=B1_L(temp_idx(1:B1_num));
            comm_rep(i).B1_beta1(1:B1_num)=B1_beta1(temp_idx(1:B1_num));
%             comm_rep(i).B1_umax(1:B1_num)=B1_umax(temp_idx(1:B1_num));
%             comm_rep(i).B1_KmFold(1:B1_num)=B1_KmFold(temp_idx(1:B1_num));
            comm_rep(i).B1_KsFold(1:B1_num)=B1_KsFold(temp_idx(1:B1_num));
        end
        if B2_num>=1
            comm_rep(i).B2_L(1:B2_num)=B2_L(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
%             comm_rep(i).B2_umax(1:B2_num)=B2_umax(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
%             comm_rep(i).B2_KsFold(1:B2_num)=B2_KsFold(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
        end
        comm_rep(i).parentnum=parentnum;
    end
    if rep_counter+comm_rep_num>comm_type_num*comm_rep_num
        comm_rep=comm_rep(1:comm_type_num*comm_rep_num-rep_counter);
    end
else
    for i=1:dil_factor
        temp_idx=find(rand_temp==i);
        B1_num=nnz(temp_idx<=B1_counter);
        B2_num=nnz(temp_idx>B1_counter);
        if B1_num>=1
            comm_rep(i).B1_L(1:B1_num)=B1_L(temp_idx(1:B1_num));
            comm_rep(i).B1_beta1(1:B1_num)=B1_beta1(temp_idx(1:B1_num));
%             comm_rep(i).B1_umax(1:B1_num)=B1_umax(temp_idx(1:B1_num));
%             comm_rep(i).B1_KmFold(1:B1_num)=B1_KmFold(temp_idx(1:B1_num));
            comm_rep(i).B1_KsFold(1:B1_num)=B1_KsFold(temp_idx(1:B1_num));
        end
        if B2_num>=1
            comm_rep(i).B2_L(1:B2_num)=B2_L(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
%             comm_rep(i).B2_umax(1:B2_num)=B2_umax(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
%             comm_rep(i).B2_KsFold(1:B2_num)=B2_KsFold(temp_idx(1+B1_num:B2_num+B1_num)-B1_counter);
        end
        comm_rep(i).parentnum=parentnum;
    end
    if rep_counter+dil_factor>comm_type_num*comm_rep_num
        comm_rep=comm_rep(1:comm_type_num*comm_rep_num-rep_counter);
    else
        comm_rep=comm_rep(1:dil_factor);
    end
end
    
    
