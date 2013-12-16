function valuefniter_case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, beta, ReturnFn, LowMemory=0,Parallel=0, Tolerance=10^(-9.0), PolIndOrVal=1,Howards=80,Verbose=0)
  #return V Policy
  
  N_d=prod(n_d)
  N_a=prod(n_a)
  N_z=prod(n_z)
  num_d=length(n_d)
  num_a=length(n_a)
  num_z=length(n_z)
  
  V0Kron=reshape(V0,N_a,N_z)

  
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
    ReturnMatrix=createreturnfnmatrix_case1(ReturnFn, n_d, n_a, n_z, d_grid, a_grid, z_grid)
    
    if Verbose==1
#      time=toc()
      println("Time to create return fn matrix: $time ")
#      tic()
    end
     
    if n_d[1]==0
      if Parallel==0
        VKron, Policy = valuefniter_case1_nod_raw(V0Kron, n_a, n_z, pi_z, beta, ReturnMatrix,Tolerance,Howards,Verbose)
#      elseif Parallel==1
        #        VKron,Policy=valuefniter_case1_nod_par1_raw(V0Kron, n_a, n_z, pi_z, beta, ReturnMatrix,Tolerance, Howards, Verbose)
      end
    else
      if Parallel==0
        VKron, Policy=valuefniter_case1_raw(V0Kron, n_d,n_a,n_z, pi_z, beta, ReturnMatrix, Tolerance, Howards, Verbose)
#      elseif Parallel==1
        #        VKron, Policy=valuefniter_case1_par1_raw(V0Kron, n_d,n_a,n_z, pi_z, beta, ReturnMatrix,Tolerance,Howards,Verbose)
      end
     end #if n_d[1]==0
    
  elseif LowMemory==1    

    if n_d[1]==0
      if Parallel==0
        #        VKron,Policy=valuefniter_case1_lowmem_nod_raw(V0Kron, n_a, n_z, a_grid, z_grid, pi_z, beta, ReturnFn, Tolerance, Howards, Verbose)
      elseif Parallel==1
        #        VKron,Policy=valuefniter_case1_lowmem_nod_par1_raw(V0Kron, n_a, n_z, a_grid, z_grid, pi_z, beta, ReturnFn,Tolerance, Howards, Verbose)
      end
    else
      if Parallel==0
        VKron, Policy=valuefniter_case1_lowmem_raw(V0Kron, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, beta, ReturnFn, Tolerance,Howards,Verbose)
      elseif Parallel==1
#        VKron, Policy=valuefniter_case1_lowmem_par1_raw(V0Kron, n_d,n_a,n_z, d_grid,a_grid,z_grid,pi_z, beta, ReturnFn,Tolerance,Howards,Verbose)
      end
    end
  end #if LowMemory==0 elseif LowMemory==1

        
    if PolIndOrVal==2
      if n_d[1]==0
        PolicyInd=copy(Policy)
        Policy=zeros(num_a,N_a,N_z)
        for a_c=1:N_a
          for z_c=1:N_z
#            temp_a=ind2grid_homemade(PolicyInd[a_c,z_c],n_a,a_grid)
            for ii=1:num_a
              Policy[ii,a_c,z_c]=temp_a[ii]
            end
          end
        end
      else #if n_d[1]==0
        PolicyInd=copy(Policy)
        Policy=zeros(num_d+num_a,N_a,N_z)
        for a_c=1:N_a
          for z_c=1:N_z
#            temp_d=ind2grid_homemade(n_d,PolicyInd[1,a_c,z_c],d_grid)
            for ii=1:num_d
              Policy[ii,a_c,z_c]=temp_d[ii]
            end
#            temp_a=ind2grid_homemade(n_a,PolicyInd[2,a_c,z_c],a_grid)
            for ii=1:num_a
              Policy[num_d+ii,a_c,z_c]=temp_a[ii]
            end
          end
        end
      end
    end #if PolIndOrVal==2



  if Verbose==1
#    time=toc()
    println("Time to solve for Value Fn and Policy: $time")
    println("Transforming Value Fn and Optimal Policy matrices back out of Kronecker Form")
#    tic()
  end
  V=reshape(VKron,n_a,n_z);
  if PolIndOrVal==1
#    Policy=UnKronPolicyIndexes_Case1(Policy, n_d, n_a, n_z)
  end
  if Verbose==1
#    time=toc()
    println("Time to create UnKron Value Fn and Policy: $time")
  end
    

  return V, Policy

end 
