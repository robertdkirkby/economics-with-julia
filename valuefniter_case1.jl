function valuefniter_case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, beta, ReturnFn, Tolerance=10^(-9.0), LowMemory=0, PolIndOrVal=1,Howards=80,Parallel=0,Verbose=0)
  #return V Policy
  
  N_d=prod(n_d);
  N_a=prod(n_a);
  N_z=prod(n_z);
  num_d=length(n_d);
  num_a=length(n_a);
  num_z=length(n_z);
  
  ##Check the sizes of some of the inputs
  # if length(n_z)==1 && n_z(1)==1
  #     if size(pi_z)~=[N_z, N_z]
  #         disp('Error: pi is not of size N_z-by-N_z')
  #     elseif size(Fmatrix)~=[n_d, n_a, n_a]
  #         disp('Error: Fmatrix is not of size [n_d, n_a, n_a, n_z]')
  #     elseif size(V0)~=[n_a]
  #         disp('Error: Starting choice for ValueFn is not of size [n_a,n_z]')
  #     end
  # else
  #     if size(pi_z)~=[N_z, N_z]
  #         disp('Error: pi is not of size N_z-by-N_z')
  #     elseif size(Fmatrix)~=[n_d, n_a, n_a, n_z]
  #         disp('Error: Fmatrix is not of size [n_d, n_a, n_a, n_z]')
  #     elseif size(V0)~=[n_a,n_z]
  #         disp('Error: Starting choice for ValueFn is not of size [n_a,n_z]')
  #     end
  # end


  if LowMemory==0
    if Verbose==1
      println("Creating return fn matrix")
#      tic()
    end
    ## CreateReturnFnMatrix_Case1_Disc creates a matrix of dimension (d and aprime)-by-a-by-z.
    # Since the return function is independent of time creating it once and
    # then using it every iteration is good for speed, but it does use a lot of memory.
    ReturnMatrix=createreturnfnmatrix_case1(ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid);
    
    if Verbose==1
#      time=toc();
      println("Time to create return fn matrix: $time ")
#      tic()
    end
    
    V0Kron=reshape(V0,N_a,N_z);
    
    #println("Starting Value Function") #No need to say this as the individual 'raw' codes will print this.
    if n_d[1]==0
      if Parallel==0
        VKron, Policy = valuefniter_case1_nod_raw(V0Kron, n_a, n_z, pi_z, beta, ReturnMatrix);
#      elseif Parallel==1
#        [VKron,Policy]=ValueFnIter_Case1_NoD_Par1_raw(Tolerance, V0Kron, n_a, n_z, pi_z, beta, ReturnMatrix, Howards, Verbose);
      end
      if PolIndOrVal==2
        PolicyInd=Policy;
        Policy=zeros(length(n_a),N_a,N_z); #NOTE: this is not actually in Kron form
        for a_c=1:N_a
          for z_c=1:N_z
#            temp_a=ind2grid_homemade(PolicyInd(a_c,z_c),n_a,a_grid);
            for ii=1:num_a
              Policy[ii,a_c,z_c]=temp_a[ii];
            end
          end
        end
      end
#    else
#      if Parallel==0
#        [VKron, Policy]=ValueFnIter_Case1_raw(Tolerance, V0Kron, n_d,n_a,n_z, pi_z, beta, ReturnMatrix,Howards,Verbose);
#      elseif Parallel==1
#        [VKron, Policy]=ValueFnIter_Case1_Par1_raw(Tolerance, V0Kron, n_d,n_a,n_z, pi_z, beta, ReturnMatrix,Howards,Verbose);
#      end
#      if PolIndOrVal==2
#        PolicyInd=Policy;
#        Policy=zeros(length(n_d)+length(n_a),N_a,N_z);
#        for a_c=1:N_a
#          for z_c=1:N_z
#            temp_d=ind2grid_homemade(n_d,PolicyInd(1,a_c,z_c),d_grid);
#            for ii=1:length(n_d)
#              Policy(ii,a_c,z_c)=temp_d(ii);
#            end
#            temp_a=ind2grid_homemade(n_a,PolicyInd(2,a_c,z_c),a_grid);
#            for ii=1:length(n_a)
#              Policy(length(n_d)+ii,a_c,z_c)=temp_a(ii);
#            end
#          end
#        end
#        clear PolicyInd
      end #if n_d[1]==0
    
#  elseif LowMemory==1    
#
#    V0Kron=reshape(V0,[N_a,N_z]);
#
#    #println("Starting Value Function") #No need to say this as the individual 'raw' codes will print this.
#    if n_d[1]==0
#      if Parallel==0
#        [VKron,Policy]=ValueFnIter_Case1_LowMem_NoD_raw(Tolerance, V0Kron, n_a, n_z, a_grid, z_grid, pi_z, beta, ReturnFn, Howards, Verbose);
#      elseif Parallel==1
#        [VKron,Policy]=ValueFnIter_Case1_LowMem_NoD_Par1_raw(Tolerance, V0Kron, n_a, n_z, a_grid, z_grid, pi_z, beta, ReturnFn, Howards, Verbose);
#      end
#      if PolIndOrVal==2
#        PolicyInd=Policy;
#        Policy=zeros(length(n_a),N_a,N_z); #NOTE: this is not actually in Kron form
#        for a_c=1:N_a
#          for z_c=1:N_z
#            temp_a=ind2grid_homemade(PolicyInd(a_c,z_c),n_a,a_grid);
#            for ii=1:length(n_a)
#              Policy(ii,a_c,z_c)=temp_a(ii);
#            end
#          end
#        end
#        clear PolicyInd
#      end
#    else
#      if Parallel==0
#        [VKron, Policy]=ValueFnIter_Case1_LowMem_raw(Tolerance, V0Kron, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, beta, ReturnFn,Howards,Verbose);
#      elseif Parallel==1
#        [VKron, Policy]=ValueFnIter_Case1_LowMem_Par1_raw(Tolerance, V0Kron, n_d,n_a,n_z, d_grid,a_grid,z_grid,pi_z, beta, ReturnFn,Howards,Verbose);
#      end
#      if PolIndOrVal==2
#        PolicyInd=Policy;
#        Policy=zeros(length(n_d)+length(n_a),N_a,N_z);
#        for a_c=1:N_a
#          for z_c=1:N_z
#            temp_d=ind2grid_homemade(n_d,PolicyInd(1,a_c,z_c),d_grid);
#            for ii=1:length(n_d)
#              Policy(ii,a_c,z_c)=temp_d(ii);
#            end
#            temp_a=ind2grid_homemade(n_a,PolicyInd(2,a_c,z_c),a_grid);
#            for ii=1:length(n_a)
#              Policy(length(n_d)+ii,a_c,z_c)=temp_a(ii);
#            end
#          end
#        end
#        clear PolicyInd
#      end
#    end
    end #if LowMemory==0 elseif LowMemory==1


  if Verbose==1
#    time=toc();
    println("Time to solve for Value Fn and Policy: $time")
    println("Transforming Value Fn and Optimal Policy matrices back out of Kronecker Form")
#    tic()
  end
  V=reshape(VKron,n_a,n_z);
  if PolIndOrVal==1
#    Policy=UnKronPolicyIndexes_Case1(Policy, n_d, n_a, n_z);
  end
  if Verbose==1
#    time=toc();
    println("Time to create UnKron Value Fn and Policy: $time")
  end
    

  return V, Policy

end 
