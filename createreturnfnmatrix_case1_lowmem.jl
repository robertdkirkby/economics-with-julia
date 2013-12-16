function createreturnfnmatrix_case1_lowmem(a,z,ReturnFn, n_d, n_a, d_gridvals, a_gridvals, Parallel=0)
  #return Fmatrix

  #Currently the 'Parallel' input is just ignored; in effect the only value it can take is 0, namely no parallelization.
  
  #If there is no d variable, just input n_d=0 and d_gridvals=0

  N_d=prod(n_d);
  N_a=prod(n_a);
  num_d=length(n_d); # In Julia you cannot current do, eg. a[1:length(a)-1], have to set l_a=length(a) and do a[1:l_a]
  num_a=length(n_a);

  if Parallel==0
  
    if N_d==0
      Fmatrix=zeros(N_a);
      for i1=1:N_a
        #aprime=a_gridvals[i1,:]; 
        Ftemp=ReturnFn(a_gridvals[i1,:],a,z);
        Fmatrix[i1]=Ftemp[1]; #Have to use Ftemp as otherwise get error "no method convert(Type{Float64}, Array{Float64,2})"
      end
    else      
      Fmatrix=zeros(N_d*N_a);
      for i1=1:N_d            
        for i2=1:N_a
          #aprime_gridvals=ind2grid_homemade(n_a,i2,a_grid);
          #i1i2=sub2ind_homemade([N_d,N_a],[i1,i2]);
          i1i2=i1+(i2-1)*N_d;
          Ftemp=ReturnFn(d_gridvals[i1,:],a_gridvals[i2,:],a,z);
          Fmatrix[i1i2]=Ftemp[1];
        end
      end
    end
    
  #elseif Parallel==1
  
  end #if Parallel==0 elseif Parallel==1


  return Fmatrix

end 