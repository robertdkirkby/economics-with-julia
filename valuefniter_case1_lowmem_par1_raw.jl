function valuefniter_case1_lowmem_raw(VKron, n_d,n_a,n_z, d_grid, a_grid, z_grid, pi_z, beta, ReturnFn, Tolerance=10^(-9.0), Howards=80,Verbose=0)
  #return VKron, Policy
  
  N_d=prod(n_d)
  N_a=prod(n_a)
  N_z=prod(n_z)
  num_d=length(n_d)
  num_a=length(n_a)
  num_z=length(n_z)

  if Verbose==1
    println("Starting Value Fn Iteration")
    tempcounter=1;
  end

  PolicyIndexes1=zeros(N_a,N_z)
  PolicyIndexes2=zeros(N_a,N_z)
  currdist=Inf

  z_gridvals=zeros(N_z,num_z);  
  for i1=1:N_z
    sub=zeros(1,num_z);
    sub[1]=rem(i1-1,n_z[1])+1;
    if num_z>1
      for ii=2:num_z-1
        sub[ii]=rem(ceil(i1/prod(n_z[1:ii-1]))-1,n_z[ii])+1;
      end
      sub[num_z]=ceil(i1/prod(n_z[1:(num_z-1)]));
      
      sub=sub+[0,cumsum(n_z[1:end-1])];
    end
    z_gridvals[i1,:]=z_grid[sub];
  end
  a_gridvals=zeros(N_a,num_a);
  for i2=1:N_a
    sub=zeros(1,num_a);
    sub[1]=rem(i2-1,n_a[1])+1;
    if num_a>1
      for ii=2:num_a-1
        sub[ii]=rem(ceil(i2/prod(n_a[1:ii-1]))-1,n_a[ii])+1;
      end
      sub[num_a]=ceil(i2/prod(n_a[1:num_a-1]));
        
      sub=sub+[0,cumsum(n_a[1:end-1])];
    end
    a_gridvals[i2,:]=a_grid[sub];
  end
  d_gridvals=zeros(N_d,num_d);
  for i3=1:N_d
    sub=zeros(1,num_d);
    sub[1]=rem(i3-1,n_d[1])+1;
    if num_d>1
      for ii=2:num_d-1
        sub[ii]=rem(ceil(i3/prod(n_d[1:ii-1]))-1,n_d[ii])+1;
      end
      sub[num_d]=ceil(i3/prod(n_d[1:num_d-1]));
    
      sub=sub+[0,cumsum(n_d[1:end-1])];
    end
    d_gridvals[i3,:]=d_grid[sub];
  end


  Ftemp_UsedForHowards=zeros(n_a,n_z);
  while currdist>Tolerance
    
    VKronold=copy(VKron)
    
    for z_c=1:N_z
      z=z_gridvals[z_c,:]
      pi_z_z=pi_z[z_c,:]
      VKron_z=zeros(N_a,1)
      PolicyIndexes1_z=zeros(N_a,1);
      PolicyIndexes2_z=zeros(N_a,1);
      Ftemp_UsedForHowards_z=zeros(n_a,1);
      #Calc the condl expectation term (except beta), which depends on z but not on control variables
#         EV_z=zeros(N_a,1); %aprime
#         for zprime_c=1:N_z
#             if pi_z(z_c,zprime_c)~=0 %multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
#                 EV_z=EV_z+VKronold(:,zprime_c)*pi_z(z_c,zprime_c);
#             end
#         end
      EV_z=VKronold.*kron(pi_z_z[1,:],ones(N_a,1));
      EV_z[isnan(EV_z)]=0; #multilications of -Inf with 0 gives NaN, this replaces them with zeros (as the zeros come from the transition probabilites)
      EV_z=sum(EV_z,2);
      entireEV_z=kron(EV_z,ones(N_d,1));
        
      for a_c=1:N_a
        #Calc the RHS
        a=a_gridvals[a_c,:]
        Fmatrix_az=createreturnfnmatrix_case1_lowmem(a,z,ReturnFn, n_d, n_a, d_gridvals, a_gridvals,0);
        entireRHS=Fmatrix_az+beta*entireEV_z; #d by aprime by 1
            
        #Calc the max and it's index
        maxindex=indmax(entireRHS);
        Vtemp=entireRHS[maxindex];
        VKron_z[a_c]=Vtemp[1];
        PolicyIndexes1_z[a_c]=rem(maxindex-1,N_d)+1;
        PolicyIndexes2_z[a_c]=ceil(maxindex/N_d);
        
        Ftemp_UsedForHowards_z[a_c,1]=Fmatrix_az[maxindex];
      end
        
      VKron[:,z_c]=VKron_z;
      PolicyIndexes1[:,z_c]=PolicyIndexes1_z;
      PolicyIndexes2[:,z_c]=PolicyIndexes2_z;
      Ftemp_UsedForHowards[:,z_c]=Ftemp_UsedForHowards_z;
    end
    
    VKrondist=VKron-VKronold; VKrondist[isnan(VKrondist)]=0;
    currdist=maximum(abs(VKrondist));
    
    if isfinite(currdist) #Use Howards Policy Fn Iteration Improvement
      for Howards_counter=1:Howards
        VKrontemp=copy(VKron)
        for z_c=1:N_z
          EVKrontemp_z=VKrontemp[PolicyIndexes2[:,z_c],:].*kron(pi_z[z_c,:],ones(N_a,1));
          EVKrontemp_z[isnan(EVKrontemp_z)]=0; #Multiplying zero (transition prob) by -Inf (value fn) gives NaN
          VKron[:,z_c]=Ftemp_UsedForHowards[:,z_c]+beta*sum(EVKrontemp_z,2);
        end
      end
    end
    
    if Verbose==1
        if rem(tempcounter,100)==0
            println("tempcounter = $tempcounter")
            println("currdist = $currdist")
        end
        tempcounter=tempcounter+1;
    end
    
  end

  return VKron, Policy

end 