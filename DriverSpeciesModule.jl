__precompile__()
module DriverSpeciesModule

using Hungarian, Missings # the hungarian algorithm used to find driver species
using MatrixNetworks
using DifferentialEquations
using ControlSystems #needed for dlqr
using StatsBase # needed for sampling in RewireNetwork

#functions that the module exports
export ER_EcologicalNetwork, DriverSpecies, SystemDynamics, System_Dynamics, f_HollingII, GetFinalState
export random_unit_vector, RewireNetwork
export GetControllerSuccess, GetControllerSucess_Rate, GetControllerSucess_Rate_OverNetworkRealizations

    """
    Returns a random ER network with:
    N = number of species,
    C = connectivity
    sigma = inter-species interaction strength

    """
    function ER_EcologicalNetwork(N::Int, C::Real, sigma:: Real)
        ## N = number of nodes
        ## C = connectivity
        ## sigma = interspecies interaction strength

        # (there are N*(N - 1) edges in a complete digraph with N nodes)
        L = floor(C* N*(N - 1)) # random edges to compute

        # initialize and fill the  AdjacencyMatrix
        A = spzeros(N,N) #initialize with zeros()

        # splitting j = 1:i-1,i+1,n avoids the if statement i==j
        for i=1:N
              for j=1:i-1
                if rand()<=C; A[i,j]=sigma*randn() end
              end
              for j=i+1:N
                if rand()<=C; A[i,j]=sigma*randn() end
              end
        end
        return(A);
    end


    """
    Returns a minimal set of driver species given the adjacency matrix A of
    the ecological network.
    A should not contain self-loops.
    Adapted from the code by Sergio Pequito
    https://la.mathworks.com/matlabcentral/fileexchange/46848-a-framework-for-structural-input-output-and-control-configuration-selection-in-large-scale-systems?s_tid=prof_contriblnk

    """
    function DriverSpecies(AdjacencyM :: SparseMatrixCSC{Float64,Int64})
        Atemp = sparse(convert(Array, AdjacencyM)) #adjacency matrix
        A = Atemp'

        # compute strongly connected components
        # note that SCC.number, SCC.map. The last one gives to which SCC belong each node
        SCC = scomponents(A)

        # Get all nodes in the i-th SCC.
        # SCCs[i] will contain all nodes in the i-th SCC.
        SCCs = Any[]; #initialize empty
        for i = 1 : SCC.number
            # find nodes in the i-th SCC, and push them to the SCCs list
            SCCs = push!(SCCs, find(SCC.map .== i))
        end

        # find the non-top likned SCCs
        nonTopLinkedSCCsIndex = Any[]; #intialize empty
        for i = 1 : SCC.number #for each SCC
            minPathToSCCj = Inf
            for j = i + 1 : SCC.number
                source = SCCs[j][1]
                target = SCCs[i][1]

                # apply breath-first search
                # d is a vector of the distance of each node to the source
                d, dt, pred = bfs(A, source)

                if d[target] .== -1
                    minPathFromiToj = Inf # d[i] = -1 denotes that the distance is inf
                else
                    minPathFromiToj = d[target]
                end

                if minPathFromiToj .< minPathToSCCj
                   minPathToSCCj = minPathFromiToj
                end
            end

            # test if there is no finite path, and in such case it is nonTopLinked
            if minPathToSCCj .== Inf
                nonTopLinkedSCCsIndex = push!(nonTopLinkedSCCsIndex,i)
            end
        end

        # Re-label of the beta (i.e.,numberNonTopLinkedSCCs) non-top linked SCCs to
        # have labels from 1 to beta.
        numberNonTopLinkedSCCs = length(nonTopLinkedSCCsIndex)
        nonTopLinkedSCCs = SCCs[nonTopLinkedSCCsIndex]

        # The bipartite graph is described by the concatenated matrix [A I]
        # where I comprises the columns with non-zero entries from the variable i to
        # x_j if x_j belongs to the non-top linked SCC i
        Imatrix = spzeros(size(A,1),numberNonTopLinkedSCCs)
        for i = 1 : numberNonTopLinkedSCCs
           Imatrix[nonTopLinkedSCCs[i], i] = 1.0
        end

        # Cost in A is inf (corresponding to "missing") in their zeros, and unitary
        # cost in there non-zero entries. Similarly, the non-zero entries of I have
        # cost 2 and zero entries cost inf (corresponding to "missing")
        # Note here we use A'
        costA = convert(Array{Union{Float64, Missings.Missing}}, abs.(sign.(A')))
        costI = convert(Array{Union{Float64, Missings.Missing}}, 2.0*abs.(sign.(Imatrix)))
        costA[find(costA .== 0.0)] = missing
        costI[find(costI .== 0.0)] = missing


        #The weighted bipartite graph is described by the concatenated matrix [costA costI]
        #weightedBipartite = convert(Array, hcat(costA, costI))
        weightedBipartite = convert(Array{Union{Float64, Missings.Missing}}, hcat(costA, costI))

        # solve the assignment problem
        # The output assignment gives a column vector of size n where the row index
        # is the state variable and the corresponding entry gives the column that as
        # been assigned. Therefore the edges in the maximum matching are given by
        # (assignment(label),label). Therefore, if an assignment(label)==0 it means
        # that there is no assignment.
        assignment, cost = hungarian(weightedBipartite)

        # The entries in assignment(label) that are greater than n (size of state
        # space) correspond to the root of an edge in E_{I,X} hence the state
        # variable indexed by the line label in assignment is a right-unmatched
        # vertex in a non-top linked SCC.
        Theta = find( assignment .> size(A,1))

        #Furthermore, the remaining right-unmatched vertices are given by
        UrMinusTheta = find( assignment .== 0.0)

        # The set of right-unmatched vertices with respect to some maximum matching of
        # B(X,X,E_{X,X}) is given by
        #Ur = hcat(Theta', UrMinusTheta')
        Ur = hcat(Theta', UrMinusTheta')

        # Gives the label of the non-top linked SCCs that have been assigned, i.e.,
        # contain a right-unmatched vertex
        nonTopLinkedSCCsAssigned = mod.(assignment[Theta], size(A,1));

        # Determine the label of a state variable for each non-top
        # linked SCC that was not assigned
        nonTopLinkedSCCsNOTassigned = setdiff(1 : numberNonTopLinkedSCCs, nonTopLinkedSCCsAssigned)

        Au = Int64[]
        for i = 1:length(nonTopLinkedSCCsNOTassigned)
            # It chooses the "first" state variable from the non-top linked SCC that
            # was not assigned. Therefore, several other alternatives are possible.
            nonTopAux = nonTopLinkedSCCs[nonTopLinkedSCCsNOTassigned[i]]
            Au = push!(Au, nonTopAux[1])
        end

        # The minimal FDIC comprising the set of right-unmatched vertices and a
        # subset of one state variable per each non-top linked SCC that was not
        # assigned.
        # mFDIC = unique([Ur Au])
        temp2 = vec(convert(Array{Int16}, Ur)) # transform Ur to a vector so it can be appended with Au
        mFDIC = unique(append!(temp2, Au))

        return(mFDIC)
    end

    """
    Rewires a network
    A = adjacency matrix of the network to rewire
    p = [0,1] probability of rewire
    """
    function RewireNetwork(A, p)
        N = size(A,1)

        row, col, s  = findnz(A) # move row, column and values  of the nonzero entries
        if p .> 0.0
                set = find(rand(nnz(A), 1) .< p) #set of indices that will be changed

                for i = 1:length(set)

                        startnode = col[set[i]];

                        index = find(col .== startnode) # find all end-nodes of start-node
                        alreadyconnected = row[index] # all nodes that start in start-node, of the form start-node -> j.
                        possibletargets = setdiff(1:N, alreadyconnected)  #nodes not already connected
                        if length(possibletargets) .>= 1
                                row[set[i]] = sample(possibletargets, 1)[1] # choose one node to replace the existing
                        end
                end
        end
        Arewired = sparse(row,col,s, N, N)

        return Arewired

    end




    """
    Returns the right-hand side of dot x = f(x) where
    x = current state
    p = vector of paramters of the form [A r theta]
    t = time (nothing to do here)

    """
    function System_Dynamics(x, p, t)
             N = length(x)

             # get (A,r, theta) from p = [A r theta]
             A = p[1:N, 1:N]
             r = p[1:N, N+1]
             theta = p[1:N, N+2]

             # if some species start with negative abundance make them zero
             pos = find(x .<= 0)
             x[pos] = 0

             # system dynamics
             dx = diagm(vec(x))*( A*f_HollingII(x, theta) + r )
             #dx =  A*x + r

             # constrian trajectories to remain positive
             dx[pos] = 0

             return(dx)
    end

    """
    Returns the functional response for Holling Type II where
            x = array(N,1) is the state
            theta  = array(N,1) is the value of the deformation
    """
    function f_HollingII(x, theta)
        x./(1.0 + theta.*x)
    end

    """
    Returns the final state of tau time units using System_Dynamics
    x0 = initial state
    tau = length of the time interval to integrate
    A = interaction matrix
    r = intrinsic growth rate vector
    theta = deformation of HollingII functional response
    """
    function GetFinalState(x0, tau, A, r, theta)
        # build the ODE problem
        #prob = ODEProblem(System_Dynamics, x0, (0.0,tau), [A r theta])
        prob = ODEProblem(System_Dynamics, x0, (0.0,tau), hcat(A, r, theta))
        # compute solution
        sol = solve(prob, reltol=1e-7, abstol=1e-7)
        #sol = solve(prob, Rodas4(autodiff=false), reltol=1e-7, abstol=1e-7)

        # find final state
        xf = sol.u[end]

        return(xf)
    end



    """
    Random unit vector in dimension d. Note that normal distribution is OK
    https://stackoverflow.com/questions/9750908/how-to-generate-a-unit-vector-pointing-in-a-random-direction-with-isotropic-dist
    """
    function random_unit_vector(d::Int)
      v = randn(d)
      return v / norm(v)
    end




   """
    Test the success of the controller with a minimal set of driver species.
    Returns a pair (ctr_sucess, m/N).
    ctr_sucess is:
     -1 if x converges to xd without control
     1 if x does not converge to xd without control, but x -> xd with control (controller succeeds)
     0 if x does not converge to xd without control, and x does not converge to xd with control (controller fails)

     m/N is: the proportion of driver species
    """
    function GetControllerSuccess(
            x0, xd,
            tau, testTime,
            MaxNumberOfSteps,
            unstableValue,
            epsilon,
            A, r, theta,
            RewiringProbability)

        N = length(x0)

        #if some species starts at zero, make it start at a positive but small value
        x0[find(x0 .<= 0)] = 0.1

        ctr_success = 0
        ## first test if the system is unstable at initial condition is unstable
        xfNC = GetFinalState(x0, testTime, A, r, theta)
        # if it is unstable, then it makes sense to apply the control
        #if (norm(xfNC - xd, Inf) .>= unstableValue) # check for instability
        # if it is unstable, then it makes sense to apply the control
        if (norm(xfNC - xd, Inf) .>= norm(x0 - xd, Inf)) # check for instability
            # find driver species
            #Atemp = copy(sparse(A))
            #Atemp[diagind(Atemp)] = 0 # Atemp = A with zero diagonal

            #Atemp_rewired = RewireNetwork(sparse(Atemp), RewiringProbability)# rewire network
            #A_rewired = sparse(Atemp_rewired - 1.0*eye(N))

            Atemp = sparse(A + 1.0*eye(N)) # A with zero diagonal
            Atemp_rewired = RewireNetwork(Atemp, RewiringProbability)
            A_rewired = Atemp_rewired - 1.0*eye(N)

            minimal_driver_species = DriverSpecies(sparse(A_rewired))
            m = length(minimal_driver_species)

            #build B matrix
            B = spzeros(N, m)
            for i = 1:m
                B[minimal_driver_species[i], i] = 1.0*rand()
            end

            # compute solution to the LQR problem
            Q = 2e4*eye(N)
            R = 1.5e-1*eye(m)
            K = dlqr(expm(A_rewired*tau), expm(A_rewired*tau)*B, Q, R)

            xprev = copy(x0)
            for k = 1:MaxNumberOfSteps # apply the iMPC
                xnext = GetFinalState(xprev - B*K*(xprev - xd), tau, A, r, theta)

                #comment this lines for silent excecution
                #print(k)
                #println(" |xnext - xd| = ", norm(xnext - xd, Inf))


                if ( norm(xnext - xd, Inf) .<= epsilon) # if control succeeded
                    ctr_success = 1
                    println("Success! |xd - x_final| = ", norm(xnext - xd, Inf))
                    break
                elseif ( norm(xnext - xd, Inf) .>= 7.0) # trajectory is too far!
                    ctr_success = 0
                    break
                else
                    pos = find(xnext .<= 0)
                    xnext[pos] = 0 # if some species is negative, make it zero
                    xprev = xnext # go for another intervention
                end
            end
            prop_DriverSpecies = m/N

        else # the system is itself stable, so we should discard this initial condition
            ctr_success = -1
            prop_DriverSpecies = -1
            println("Stable by itself...")
        end
        return ctr_success, prop_DriverSpecies
    end


    """
    Test the sucess rate of the controller with a minimal set of driver species
    over a collection of initial conditions (parametrized by distance).
    Returns a pair (sucess_rate, mean proportion of driver species)
    """
    function GetControllerSucess_Rate(
            xd,
            distance,
            NumberOfInitialConditions,
            tau, testTime,
            MaxNumberOfSteps,
            unstableValue,
            A, r, theta,
            RewiringProbability
            )

        N = length(xd)

        # two vectors to accumulate the result of GetControllerSucess
        ctrl_success = Int64[]
        prop_DriverSpecies = Float64[]

        for k =1:NumberOfInitialConditions
            x0 = xd + distance*random_unit_vector(N) # generate initial condition
            epsilon = 0.1*norm(x0 - xd, Inf)
            sucess, proportion_DriverSpecies = GetControllerSuccess(x0, xd, tau, testTime, MaxNumberOfSteps, unstableValue, epsilon,  A, r, theta, RewiringProbability)
            ctrl_success = push!(ctrl_success, sucess)
            prop_DriverSpecies = push!(prop_DriverSpecies, proportion_DriverSpecies)
        end

        # return the mean of initial conditions that were not stable by themselves
        return mean(ctrl_success[find(ctrl_success .>= 0)]), mean(prop_DriverSpecies[find(prop_DriverSpecies .>= 0)])
    end


    """
    Test the sucess rate of the controller with a minimal set of driver species,
    over a collection of initial conditions (parametrized by distance) and
    Number of Networks (parametrized by N, C, sigma).
    Returns a pair of vectors (success_rate_vector, mean_prop_driverspecies_vector)
    """
    function GetControllerSucess_Rate_OverNetworkRealizations(
        distance,
        NumberOfInitialConditions,
        NumberOfNetworks,
        N, C, sigma,
        RewiringProbability)

        delta = 1.0 # inter-species interactions strength
        theta_max = 0.05 # maximum value for the deformation
        tau = 0.1 #time between interventions


        testTime = 1 #time interval for integration to test if the dynamcis is unstable
        unstableValue = 1 # deviaton from the desired state that will be considered as evidence of instability
        MaxNumberOfSteps = 700 # maximum number of interventions for the MPC control

        xd = ones(N,1) # desired state

        success_rate_vector = Float64[]
        mean_prop_driverspecies_vector = Float64[]
        #@time
        for k = 1:NumberOfNetworks
            ## Generate the parameters (A, r, theta)
            ##
            # adjacency matrix of the ecological network, without self-loops
            A0 = ER_EcologicalNetwork(N, C, sigma)
            A = convert(Array{Float64}, A0 - delta*eye(N)) # the A matrix is obtained by add diagonal terms
            # Set value of the deformation theta
            theta = theta_max*rand(N,1)
            #to compute the growth rate, we select it such that xd = ones(N,1) is an equilibrium
            r = -A*f_HollingII(xd,theta)

            ## Compute sucess rate
            ##
            success_rate, mean_prop_driverspecies =
            GetControllerSucess_Rate(xd, distance, NumberOfInitialConditions, tau, testTime, MaxNumberOfSteps, unstableValue, A, r, theta, RewiringProbability)

            success_rate_vector = push!(success_rate_vector, success_rate)
            mean_prop_driverspecies_vector = push!(mean_prop_driverspecies_vector, mean_prop_driverspecies)
        end

        return success_rate_vector, mean_prop_driverspecies_vector
    end

end
