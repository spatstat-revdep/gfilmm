using LinearAlgebra
using Distributions # to sample chi²

## similar to the Matlab orth() function
## -->> faire qr(A).Q
function orth(A)
  if (isempty(A))
    retval = []
  else
    (U, S, V) = svd(A)
    (rows, cols) = size(A)
    tol = maximum(size(A)) * S[1] * eps()
    r = sum(S .> tol)
    if (r > 0)
      retval = -U[:, 1:r]
    else
      retval = zeros(rows, 0)
    end
  end
  return (retval)
end
## sort array M with respect to colon j
sortbycol = function (M, j)
  local col, perm
  col = M[:, j]
  perm = sortperm(col)#[2]
  return (M[perm, :])
end

function fid_vertex0(
  VT1::Array{R},
  CC1::Array{Int64},
  VTsum::Vector{R},
  U::R,
  L::R,
  m::Int64,
  dim::Int64,
  k::Int64,
  n::Int64,
) where {R<:Real}
  # VT1 : matric 4x16
  # CC1: idem integer entries
  # VTsum: colonne 16
  # U,L : nbm
  # m,dim,k,n : integers
  m = length(VTsum) # ?
  whichVT = zeros(Bool, m, 2)
  whichVT[:, 1] = VTsum .>= L
  whichVT[:, 2] = VTsum .<= U
  whichl = findall(whichVT[:, 1])  #vertices that satisfy lower constraint
  whichu = findall(whichVT[:, 2])  #vertices that satisfy upper constraint
  both = findall(all(whichVT, dims=2)[:, 1]) ## -->> findall(all(whichVT, dims=2))
  ##       -> cartesian index
  ##      findall(all(whichVT, dims=2)[:,1])
  checkl = findall(.!whichVT[:, 1])
  checku = findall(.!whichVT[:, 2])
  CCtemp = zeros(dim, 0) # ou utiliser peut-?tre cell()
  VTtemp = zeros(dim, 0)
  vert = 0
  CA = CC1[:, checkl]
  CB = CC1[:, whichl]
  if !isempty(checkl) #i.e. they do not all satisfy the lower constraint
    INT = zeros(Int64, 2 * n, length(checkl))
    for ll = 1:length(checkl) #check lower constraints first
      INT[CA[:, ll], ll] .= 1
    end
    for ii = 1:length(whichl)
      INT2 = INT[CB[:, ii], :]
      #print(INT2)
      use = findall(vec(sum(INT2, dims=1) .== dim - 1)) # sum sur colonne
      #print(use)
      for dd = 1:length(use)
        inter = CB[findall(INT2[:, use[dd]] .== 1), ii]  #this will be intersection
        #print(inter)
        vert = vert + 1
        CCtemp = [CCtemp [inter; k + n]]  # ATTENTION A VERIFIER need to add n indicating lower constraint
        lambda = (L - VTsum[whichl[ii]]) /
          (VTsum[checkl[use[dd]]] - VTsum[whichl[ii]])
        VTtemp =
          [VTtemp lambda * VT1[:, checkl[use[dd]]] + (1.0 - lambda) * VT1[:, whichl[ii]]] # il faut coller en colonnes
      end
    end
  end
  CA = CC1[:, checku]
  CB = CC1[:, whichu]
  if !isempty(checku) #i.e. they do not all satisfy the lower constraint
    INT = zeros(2 * n, length(checku))
    for ll = 1:length(checku) #check lower constraints first
      INT[CA[:, ll], ll] .= 1
    end
    for ii in 1:length(whichu)
      #print("\nINT2 below\n")
      INT2 = INT[CB[:, ii], :]
      use = findall((sum(INT2, dims=1) .== dim - 1)[1, :])
      #print("\n************use*********:\n")
      #print(use)
      #print("\n")
      #print(CB)
      for dd in 1:length(use)
        #print("\naaaaaaa\n")
        #print(use[dd])
        #print("\n")
        #print(findall((INT2[:, use[dd]] .== 1)[:,1]))
        #print("\nxxx\n")
        inter = CB[findall((INT2[:, use[dd]] .== 1)[:,1]), ii]
        #print(inter)  #this will be intersection
        vert = vert + 1
        CCtemp = [CCtemp [inter; k]]  #need to add n indicating lower constraint
        lambda = (U - VTsum[whichu[ii]]) / (VTsum[checku[use[dd]]] - VTsum[whichu[ii]])
        #print("lambda\n")
        VTtemp =
          [VTtemp lambda * VT1[:, checku[use[dd]]] + (1.0 - lambda) * VT1[:, whichu[ii]]]
        #print("VTtemp\n")
      end
    end
  end
  if !isempty(both)
    for ll in both # both' ??
      vert = vert + 1
      CCtemp = [CCtemp CC1[:, ll]]
      VTtemp = [VTtemp VT1[:, ll]]
    end
  end
  return (VTtemp, CCtemp, vert)
end

function fid_vertex(
  VTi::Array{R},
  CCi::Array{Int64},
  VTsum::Vector{R},
  Uk::R,
  Lk::R,
  p,
  dim::Int64,
  k::Int64,
  n::Int64,
) where {R<:Real}
  p = length(VTsum)
  lwr = VTsum .>= Lk
  upr = VTsum .<= Uk
  whichl = findall(lwr)  #vertices that satisfy lower constraint
  whichu = findall(upr)  #vertices that satisfy upper constraint
  both = findall(lwr .& upr)
  checkl = findall(.!lwr)
  checku = findall(.!upr)
  lcheckl = length(checkl)
  lchecku = length(checku)
  lwhichl = p - lcheckl
  lwhichu = p - lchecku
  #
  CCtemp = zeros(Int64, dim, 0)
  VTtemp = zeros(Float64, dim, 0)
  vert = 0
  #
  CA = CCi[:, checkl]
  CB = CCi[:, whichl]
  if lcheckl != 0
    INT = falses(2 * n, lcheckl)
    INT[CA, :] .= 1 #check lower constraints first
    #
    VTsum_cl = VTsum[checkl]
    VT1_cl = VTi[:, checkl]
    VTsum_wl = VTsum[whichl]
    VT1_wl = VTi[:, whichl]
    #
    for ii in 1:lwhichl
      INT2 = INT[CB[:, ii], :]
      use = findall(count(==(1), INT2, dims=1)[1, :] .== dim-1)
      #print(use)
      #print("\n")
      vert = vert + length(use)
      for dd in use
        inter = CB[findall(INT2[:, dd] .== 1), ii]  #this will be intersection
        CCtemp = hcat(CCtemp, [inter; k+n])  # ATTENTION A VERIFIER need to add n indicating lower constraint
        lambda = (Lk - VTsum[whichl[ii]]) / (VTsum[checkl[dd]] - VTsum[whichl[ii]])
        VTtemp = hcat(
          VTtemp,
          lambda * VTi[:, checkl[dd]] + (1.0 - lambda) * VTi[:, whichl[ii]]
        )
      end
    end
  end
  CA = CCi[:, checku]
  CB = CCi[:, whichu]
  if lchecku != 0
    INT = falses(2 * n, lchecku)
    INT[CA, :] .= 1 #check upper constraints
    #
    VTsum_cu = VTsum[checku]
    VT1_cu = VTi[:, checku]
    VTsum_wu = VTsum[whichu]
    VT1_wu = VTi[:, whichu]
    #
    for ii in 1:lwhichu
      INT2 = INT[CB[:, ii], :]
      use = findall(count(==(1), INT2, dims=1)[1, :] .== dim-1)
      vert = vert + length(use)
      for dd in use
        inter = CB[findall(INT2[:, dd] .== 1), ii]  #this will be intersection
        CCtemp = hcat(CCtemp, [inter; k])  # ATTENTION A VERIFIER need to add n indicating lower constraint
        lambda = (Uk - VTsum[whichu[ii]]) / (VTsum[checku[dd]] - VTsum[whichu[ii]])
        VTtemp = hcat(
          VTtemp,
          lambda * VTi[:, checku[dd]] + (1.0 - lambda) * VTi[:, whichu[ii]]
        )
      end
    end
  end
  if !isempty(both)
    CCtemp = hcat(CCtemp, CCi[:, both])
    VTtemp = hcat(VTtemp, VTi[:, both])
    vert = vert + length(both)
  end

  return (VTtemp, CCtemp, vert)
end

function approx(x::R, n::Int64) where {R<:Real}
  return round(x * 10^n) / 10^n
end

function fid_sample(VT2::Vector{R}, VTsum::Vector{R}, U::R, L::R) where {R<:Real}
  high = findall(VT2 .> 0.0)
  low = findall(VT2 .< 0.0)
  zero = findall(VT2 .== 0.0)
  if (!isempty(low) && !isempty(high)) || !isempty(zero) #some positive and negative sigs
    UU = map(sign, U .- VTsum)
    LL = map(sign, L .- VTsum)
    SS = map(sign, VT2)
    whichnot = findall(VT2 .!= 0)
    which = findall(VT2 .== 0)
    if length(zero) == length(VT2) #all are zero
      if !all(VTsum .>= U) && !all(VTsum .<= L) #do not satisfy constaints
        #changed to be strict inequalities 2/24 due to precision of
        #Matlab
        MAX = Inf
        MIN = -Inf
        temp = 0  #let Z's be anything, but give weight = 0
      else #satisfy constraints...Z can be anything
        MAX = Inf
        MIN = -Inf
        temp = 1
      end
    elseif (
      !any(SS .== -1) && !any(UU[which] .== -1) && !any(LL[which] .== -1)
    ) || (!any(SS .== 1) && !any(UU[which] .== 1) && !any(LL[which] .== 1))
      Us = (U .- VTsum[whichnot]) ./ VT2[whichnot]
      Ls = (L .- VTsum[whichnot]) ./ VT2[whichnot]
      MAX = Inf
      MIN = min(minimum(Us), minimum(Ls))
      temp = 1 - (atan(MIN) / pi + 0.5)  #weight adjustment
    elseif (
      !any(SS .== 1) && !any(UU[which] .== -1) && !any(LL[which] .== -1)
    ) || (
      !any(SS .== -1) && !any(UU[which] .== 1) && !any(LL[which] .== 1)
    )
      Us = (U .- VTsum[whichnot]) ./ VT2[whichnot]
      Ls = (L .- VTsum[whichnot]) ./ VT2[whichnot]
      MAX = max(maximum(Us), maximum(Ls))
      MIN = -Inf
      temp = atan(MAX) / pi + 0.5  #weight adjustment
    else
      if isempty(high)
        Hmax = -Inf
        Hmin = Inf
      else
        HUs = (U .- VTsum[high]) ./ VT2[high]
        HLs = (L .- VTsum[high]) ./ VT2[high]
        Hmax = max(maximum(HUs), maximum(HLs))
        Hmin = min(minimum(HUs), minimum(HLs))
      end
      if isempty(low)
        Lmax = -Inf
        Lmin = Inf
      else
        LUs = (U .- VTsum[low]) ./ VT2[low]
        LLs = (L .- VTsum[low]) ./ VT2[low]
        Lmax = max(maximum(LUs), maximum(LLs))
        Lmin = min(minimum(LUs), minimum(LLs))
      end
      if 0 <= approx(Lmin - Hmax, 12)    ## roundn indisponible
        bottom_pos = -Inf
        top_pos = Hmax
        bottom_neg = Lmin
        top_neg = Inf
      elseif approx(Hmin - Lmax, 12) >= 0
        bottom_pos = Hmin
        top_pos = Inf
        bottom_neg = -Inf
        top_neg = Lmax
      else
        bottom_pos = -Inf
        top_pos = Inf
        bottom_neg = -Inf
        top_neg = Inf
      end
      #######################################################################
      if top_pos == Inf
        Pprob = 1.0 - (atan(bottom_pos) / pi + 0.5)
      elseif bottom_pos == -Inf
        Pprob = atan(top_pos) / pi + 0.5
      end
      if top_neg == Inf
        Nprob = 1.0 - (atan(bottom_neg) / pi + 0.5)
      elseif bottom_neg == -Inf
        Nprob = atan(top_neg) / pi + 0.5
      end
      temp = Pprob + Nprob
      Pprob = Pprob / temp
      Nprob = Nprob / temp
      if rand() <= Nprob #choose negative range
        MIN = bottom_neg
        MAX = top_neg
      else #choose positive range
        MIN = bottom_pos
        MAX = top_pos
      end
    end
    #######################################################################
    y = atan(MAX) / pi + 0.5  #cdf
    x = atan(MIN) / pi + 0.5
    u = x + (y - x) * rand()
    ZZ = tan(pi * (u - 0.5))  #Inverse cdf
    wt = exp(-ZZ^2 / 2.0) * (1.0 + ZZ^2) * temp
  else
    Us = (U .- VTsum) ./ VT2
    Ls = (L .- VTsum) ./ VT2
    MAX = max(maximum(Us), maximum(Ls))
    MIN = min(minimum(Us), minimum(Ls))
    y = atan(MAX) / pi + 0.5 #cdf
    x = atan(MIN) / pi + 0.5
    u = x + (y - x) * rand()
    ZZ = tan(pi * (u - 0.5)) #Inverse cdf
    wt = exp(-ZZ^2 / 2.0) * (1.0 + ZZ^2) * (y - x)
  end
  return (ZZ, wt, MIN, MAX)
end


#srand(421)
gfilmm = function (data::Array{R,2}, FE::Array{R,2}, RE::Array{Int64,2}, N::Int64, thresh::R) where {R<:Real}
  #global AAA, AAAA, KKK, OOO2, COOO2, ZZZ1, MMMAT, ESSS, ZZZZZ, iiii, ZZZZ
  #VERTEX = []
  #WEIGHT = []
  #ESS = []
  #p = []
  #r = []
  local WT, VERTEX, WEIGHT, VTtemp, VCtemp

  n = size(data, 1)
  random_design = 1
  Y = data
  if size(data, 2) == 2
    L = Y[:, 1]
    U = Y[:, 2]
  else
    error("ERROR:  Check the form of data:  should be n-by-2")
  end
  if count(U .< L) > 0
    error("ERROR:  Check the form of data:  data(i 1)<data(i 2)  for i = 1:n")
  end
  break_point = 10  # length of steps before recycling
  if !isempty(FE)
    fe = size(FE, 2)
  else
    fe = 0  #only error component
  end

  ########--------SET-UP RANDOM EFFECTS
  if random_design == 1       # RE is declared.
    ree = size(RE, 2) + 1
    E = zeros(Int64, ree)
    E[ree] = n
    for i in 1:(ree-1)
      E[i] = length(unique(RE[:, i]))
    end
    ESUM = cumsum(E)
    RE2 = hcat(RE, 1:n)   #Adds the error effect
    RE = zeros(Int64, n, 0)
    for i in 1:ree #Builds an indicator matrix for the effects
      re_levels = unique(RE2[:, i])
      for j in 1:E[i]
        temp1 = RE2[:, i] .== re_levels[j]
        temp2 = zeros(Int64, n)
        temp2[temp1] .= 1
        RE = hcat(RE, temp2)
      end
    end
  elseif random_design == 2       #V and Levels are declared
    RE = [V eye(n)]
    E = [Levels n]
    ESUM = cumsum(E)
    ree = length(E)
    RE2 = zeros(n, ree)  # assigns the level per effect of each observation
    for i = 1:n
      for j = 1:ree
        temp = find(RE[i, sum(E[1, 1:j])-E[1, j]+1:sum(E[1, 1:j])] .== 1)
        if isempty(temp) == 0
          RE2[i, j] = temp
        end
      end
    end
  else
    RE = eye(n)
    RE2 = [1:n]'
    E = n
    ESUM = n
    ree = 1
  end

  dim = fe + ree  #dimension of the space
  rep = 1

  ########--------SET-UP ALGORITHM OBJECTS
  Z = Vector{Array{Float64,2}}(undef, ree)     #Particles
  weight = Vector{Array{Float64,2}}(undef, ree)   #Weights
  ESS = N * ones(Float64, n)     #Effective sample size
  VC = zeros(Int64, N)       #Number of vertices
  VT = Vector{Array{Float64,2}}(undef, N)        #Verticies
  CC = Vector{Array{Int64,2}}(undef, N)        #Constraints

  ########--------SAMPLE ALL Z's/SET-UP WEIGHTS

  A = Vector{Array{Float64,2}}(undef, N)
  for j in 1:N
    A[j] = zeros(Float64, n, 0)
  end

  for i in 1:ree  #Sample all the Z's and set-up the weights
    Z[i] = randn(E[i], N)
    weight[i] = ones(Float64, E[i], N)
    for j in 1:N
      A[j] = hcat(A[j], RE[:, ESUM[i]-E[i]+1:ESUM[i]] * Z[i][:, j])
    end
  end

  C = zeros(Int64, 0) #Initial constraints selected
  AA = hcat(vcat(FE, -FE), vcat(A[1], -A[1]))
  AT = zeros(Float64, 0, size(AA, 2))
  r = 0
  if rank(AA) != dim
    error("Design is not of full rank")
  else
    for k in 1:n
      A_temp = vcat(AT, AA[k, :]')
      if rank(A_temp) > r
        AT = A_temp
        r += 1
        C = vcat(C, k)
      end
    end
  end

  K = setdiff(1:n, C)
  K_start = C
  Z[ree][K, :] .= 0  #remove error particles not selected in intitialization
  # K will be the values left to be sampled

  ########--------FIND INITIAL VERTICES

  USE = repeat(C', outer = (2^dim, 1))  #all combinations of I
  for j in 1:dim
    if j == 1
      USE[(2^(dim-1)+1):(2^dim), 1] = USE[1:(2^(dim-1)), 1] .+ n
    else
      temp = 2^(j - 2) + 1
      while temp <= 2^dim
        for ii in 1:(2^(j-2))
          USE[temp, j] += n
          temp += 1
        end
        temp = temp + 2^(j - 2)
      end
    end
  end

  for i in 1:N
    V = zeros(Float64, dim, 2^dim)
    temp = 0   # ???
    for ii in 1:(2^dim) #number of vertices
      II = USE[ii, :]  #constraints to use are by row
      AA = hcat(vcat(FE, -FE), vcat(A[i], -A[i])) # ? sortir de cette boucle
      b = vcat(U, -L) # ça aussi
      AA = AA[II, :]
      b = b[II]
      V[:, ii] = AA \ b
    end
    VT[i] = V
  end
  for k in 1:N
    CC[k] = USE' #reshape(USE', dim, 2^dim)  #constraints are the same for all N particles
  end
  VC = 2^dim * ones(Int64, N)

  ########--------MAIN ALGORITHM
  K_n = convert(Int64, ceil(length(K) / break_point))
  K_temp = Vector{Vector{Int64}}(undef, K_n)
  for i in 1:(K_n-1)
    K_temp[i] = K[((i-1)*break_point+1):(i*break_point)]
  end
  K_temp[K_n] = K[((K_n-1)*break_point+1):end]

  K1 = Vector{Int64}(undef, 0)
  for k_n in 1:K_n
    K1 = vcat(K1, K_temp[k_n])
    for k in K_temp[k_n]
      effect = ree  #only have the error term left
      level = k  #the level is the observation we are on
      if k_n > 1
        for i in 1:ree
          Z[i] = vcat(Z[i], randn(E[i] - size(Z[i], 1), N))
        end
      end
      for i in 1:N
        #m = convert(Int, VC[1, i])
        m = VC[i]
        VT1 = VT[i]  #vertex values
        VT2 = VT1[fe+effect, :]  #value to be resampled
        # remove value to be resampled
        keepcols = setdiff(1:size(VT1, 1), fe + effect)
        VT1 = VT1[keepcols, :]  #
        if fe > 0
          Z1 = FE[k, :]  #Assign fixed effect values
        else
          Z1 = zeros(Float64, 0)
        end
        for j in 1:ree
          if RE2[k, j] == 0 # -->> it's never 0, I think
            Z1 = [Z1 0]
          else
            Z1 = vcat(Z1, Z[j][RE2[k, j], i])
          end
        end
        Z1 = Z1[keepcols]   #remove column of effect to be sampled
        VTsum = (Z1' * VT1)'
        (ZZ, wt) = fid_sample(VT2, VTsum, U[k], L[k]) ###Sample -
        if wt == 0.0
        #  error("wt 0")
        end
        Z[effect][k, i] = ZZ
        weight[effect][k, i] = wt # -->> only one 'weight' object is used !
        VTsum += Z[effect][k, i] * VT2
        VT1 = VT[i]
        CC1 = CC[i]
        #        (VTtemp, CCtemp, vert) = fid_vertex(VT1, CC1, VTsum, U[k], L[k], dim, k, n) # -->> I removed the 'm'
        #        VC[i] = vert
        #        CC[i] = CCtemp
        #        VT[i] = VTtemp
        (VT[i], CC[i], VC[i]) = fid_vertex0(VT1, CC1, VTsum, U[k], L[k], m, dim, k, n) # -->> I removed the 'm'
        VTtemp = VT[i]
        VCtemp = VC[i]
        if VC[i] == 0 #check
          print("\nVC[i] zero\n")
          # println(CC)
          # println(VT1)
          # println(CC1)
          # println(VTsum)
          # println([U[k] L[k]])
          # println(i)
          # println(wt)
          #error("STOP")
          weight[effect][k, i] = 0
        end
      end
      # ZZZZ[k] = Z[effect, 1] -->> not used
      WT = vec(prod(weight[ree], dims = 1)) #only last re is restricted
      # KKK = K_n - k_n -->> not used
      if sum(WT) == 0#
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("!!!   Error: possible underflow. !!!")
        print("!!!   Potential resolutions:      !!!")
        print("!!!   1.  Widen data intervals   !!!")
        print("!!!   2.  Center data            !!!")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return ()
      end
      WT /= sum(WT)
      ESS[k] = 1 / (WT' * WT)
      # ESSS[k] = ESS[k, 1] -->> notused
      #---------------------------------------------------------------Resample
      if ESS[k] < thresh && k < K[end] #- 1 # -->> pourquoi -1 ?
        u = zeros(Float64, N)
        N_sons = zeros(Int64, N)
        # generate the cumulative distribution
        dist = N * cumsum(WT)
        aux = rand()    # sample uniform rv in [0 1]
        u = aux .+ (0:(N-1))
        #u = u ./ N -->> I multiplied 'dist' by N instead
        j = 1
        for i in 1:N
          while (u[i] > dist[j])
            j += 1
          end
          N_sons[j] += 1
        end
        #KEEP[k,:]=N_sons
        #        II = setdiff(1:n, union(1:k, C))
        #        JJ = setdiff(1:n, II)
        JJ = union(1:k, K_start)
        ZZ = Vector{Array{Float64,2}}(undef, ree)
        for ii in 1:ree
          ZZ[ii] = zeros(Float64, size(Z[ii], 1), 0) # initialisation matrice vide
        end
        VCVC = zeros(Int64, N)
        CCCC = Vector{Array{Int64,2}}(undef, N)
        VTVT = Vector{Array{Float64,2}}(undef, N)
        for i in 1:N
          if N_sons[i] > 0
            VCtemp = VC[i] * ones(Int64, 1, N_sons[i])
            Ztemp = Vector{Array{Float64,2}}(undef, ree)
            VTtemp = Vector{Array{Float64,2}}(undef, N_sons[i])
            copy = N_sons[i] - 1  # # to be resampled
            for ii in 1:N_sons[i]
              VTtemp[ii] = VT[i]  #copy original vertices
            end
            for ii in 1:ree
              Ztemp[ii] = repeat(Z[ii][:, i], outer = (1, N_sons[i]))  #copy Z
            end
            ##########################################################
            if copy > 0
              for rr in 1:rep # -->> ? rep = 1 (declared at the beginning)
                ord = sortbycol(hcat(1:ree, rand(ree)), 2)
                ord = ord[:, 1]'  #Order to resample.  Each re will be resampled.
                ord = convert(Vector{Int64}, vec(ord))
                for kk in ord
                  for ii in 1:copy
                    CO = RE[JJ, :]
                    XX = zeros(Float64, n, 0)
                    for jj in 1:ree
                      XX = hcat(XX, RE[:, (ESUM[jj]-E[jj]+1):ESUM[jj]] * Ztemp[jj][:, ii])
                    end
                    XX = XX[JJ, :]
                    XX = XX[:, setdiff(1:size(XX, 2), kk)] #remove column of effect to be resampled
                    temp = RE2[JJ, kk] .!= 0  # -->> USELESS ? find which levels of kk have been sampled
                    Z1 = Ztemp[kk][unique(RE2[JJ[temp], kk]), ii]  #Z being resampled
                    CO2 = RE[JJ, (ESUM[kk]-E[kk]+1):ESUM[kk]]
                    level0 = findall(vec(all(CO2 .== 0, dims=1))) # findall(vec(sum(CO2, dims = 1)) .== 0) #levels not sampled yet # -->> TO REPLACE WITH count
                    Z0 = findall(Z1 .== 0.0)  #Z's not sampled
                    Z00 = setdiff(1:length(Z1), Z0) #These are the levels with Z for effect kk
                    CO2 = CO2[:, setdiff(1:size(CO2, 2), level0)]
                    Z1 = Z1[Z00]#setdiff(Z1, Z1[Z0]) # -->> was setminus(Z0, Z1) --- maybe use setdiff!
                    if fe > 0
                      XX = hcat(FE[JJ, :], XX)
                    end
                    MAT = hcat(-XX, CO2)
                    if rank(MAT) < size(MAT, 2) # -->> looks like it is 1, not 2
                      NUL = nullspace(MAT) #
                      n1 = NUL[1:size(NUL, 1) - size(CO2, 2), :]
                      n2 = NUL[(size(NUL, 1) - size(CO2, 2)+1):end, :]
                      #O2 = qr(n2).Q[1:size(n2, 1), 1:size(n2, 2)]
                      O2 = orth(n2)
                      B = CO2 * O2
                      O1 = XX \ B
                      a = O2' * Z1
                      b = sqrt((Z1 - O2 * a)' * (Z1 - O2 * a))
                      tau = (Z1 - O2 * a) / b
                      #AAA = length(Z1) - rank(O2)
                      bb = rand(Chisq(length(Z1) - rank(O2))) # -->> rank O2 is ncol or nrow
                      #                      bb = sqrt(randchi2(max(size(Z1)) - rank(O2)))
                      bbb = b/bb #
                      aa = randn(rank(O2), 1) # -->> idem rank O2
                      Ztemp[kk][Z00, ii] = O2 * aa + bb * tau
                      vert = setdiff(1:dim, fe + kk)
                      for jj in 1:VC[i]
                        check1 =
                          XX * VTtemp[ii][vert, jj] + VTtemp[ii][fe+kk, jj] * CO2 * Z1
                        VTtemp[ii][vert, jj] -=
                          VTtemp[ii][fe+kk, jj] * O1 * (bbb * aa - a)
                        VTtemp[ii][fe+kk, jj] *= bbb
                        check2 = XX * VTtemp[ii][vert, jj] +
                          VTtemp[ii][fe+kk, jj] * CO2 * (O2 * aa + bb * tau)
                      end
                    else
                      b = sqrt(Z1' * Z1)
                      tau = Z1 / b
                      bb = rand(Chisq(length(Z1)))
                      Ztemp[kk][Z00, ii] = bb * tau
                      vert = setdiff(1:dim, fe + kk)
                      for jj in 1:VC[i]
                        VTtemp[ii][fe+kk, jj] *= b / bb
                      end
                    end
                  end
                end
              end
            end
            for ii in 1:ree
              ZZ[ii] = hcat(ZZ[ii], Ztemp[ii])  # ? c'est pareil que ZZ[ii,1]=2temp[ii,1]
            end
            VCVC[(sum(N_sons[1:(i-1)])+1):sum(N_sons[1:i])] = VCtemp
            d = sum(N_sons[1:(i-1)])
            for kk in 1:N_sons[i]
              VTVT[kk+d] = VTtemp[kk]
              CCCC[kk+d] = CC[i]
            end
          end
        end
        Z = ZZ
        VT = VTVT
        VC = VCVC
        CC = CCCC
        weight[ree] = ones(Float64, E[ree], N)  #assign weights of error matrix to 1
      end
    end  #ends reampling for k=K1
    #----------------------------------------------------determine signs
    signs = zeros(Int64, ree, N)
    for i in 1:N
      for j in 1:ree
        positive = VT[i][fe+j, :] .> 0
        negative = VT[i][fe+j, :] .< 0
        if count(positive) == VC[i] #i.e. all are positive
          signs[j, i] = 1
        elseif count(negative) == VC[i] #i.e. all are negative
          signs[j, i] = -1
        end
      end
    end
    #----------------------------------------------------FINAL RESAMPLE
    ZZ = Vector{Array{Float64,2}}(undef, ree)
    VTVT = Vector{Array{Float64,2}}(undef, N)
    #ZZZZZ = Z -->> not used

    n1 = sort(vcat(K1, K_start))
    nn = Vector{Vector{Int64}}(undef, ree)
    for ii in 1:ree
      nn[ii] = unique(RE2[n1, ii])
      ZZ[ii] = zeros(Float64, length(nn[ii]), 0)
    end

    for i in 1:N
      # iiii = i -->> not used
      Ztemp = Vector{Vector{Float64}}(undef, ree)
      VTtemp = VT[i]
      for ii in 1:ree
        #nn[ii] = unique(RE2[n1, ii])
        Ztemp[ii] = Z[ii][nn[ii], i]  #copy Z
        #if true #i == 1 # initialisation que pour i=1
        #  ZZ[ii] = zeros(Float64, length(Ztemp[ii]), 0)
        #end
      end
      ord = sortbycol(hcat(1:ree, rand(ree)), 2)
      ord = ord[:, 1]'  #Order to resample.  Each re will be resampled.
      ord = convert(Array{Int,2}, ord)
      for kk in ord
        CO = RE
        XX = zeros(Float64, size(RE, 1), 0)
        eff = 1:ree
        eff = setdiff(eff, kk)
        for jj in eff
          XX = hcat(XX, RE[:, (ESUM[jj]-E[jj]+1):(ESUM[jj]-E[jj]+length(nn[jj]))]) *
              Ztemp[jj]
        end
        #temp=find(RE2[:,kk].!=0)  #find which levels of kk have been sampled
        Z1 = Ztemp[kk]  #Z being resampled
        CO2 = RE[:, (ESUM[kk]-E[kk]+1):(ESUM[kk]-E[kk]+length(nn[kk]))]  # c'?tait E[:,kk] mais une seule ligne
        level0 = findall(vec(all(CO2 .== 0, dims=1))) #levels not sampled yet
        Z0 = findall(Z1 .== 0)  #Z's not sampled
        Z00 = setdiff(1:length(Z1), Z0) #These are the levels with Z for effect kk
        CO2 = CO2[:, setdiff(1:size(CO2, 2), level0)]
        Z1 = Z1[Z00] # -->> same as setdiff(Z1, Z1[Z0])
        if fe > 0
          XX = hcat(FE, XX)
        end
        MAT = hcat(-XX, CO2)
        #MMMAT = MAT
        #COOO2 = CO2
        if rank(MAT) < size(MAT, 2) # -->> was 2
          NUL = nullspace(MAT) #
          n1 = NUL[1:(size(NUL, 1) - size(CO2, 2)), :]
          n2 = NUL[(size(NUL, 1) - size(CO2, 2) + 1):end, :]
          O2 = orth(n2)
#          O2 = qr(n2).Q[1:size(n2, 1), 1:size(n2, 2)]  # rank(O2) = ncol(NUL) ?
          B = CO2 * O2
          O1 = XX \ B
          a = O2' * Z1
          b = sqrt((Z1 - O2 * a)' * (Z1 - O2 * a))
          tau = (Z1 - O2 * a) / b
          #OOO2 = O2
          #ZZZ1 = Z1
          #AAAA = max(size(Z1)) - rank(O2)
          bb = rand(Chisq(length(Z1) - rank(O2)))
          bbb = b/bb #
          aa = randn(rank(O2), 1)
          Ztemp[kk][Z00] = O2 * aa + bb * tau
          vert = setdiff(1:dim, fe + kk)
          for jj = 1:VC[i]
            VTtemp[vert, jj] -= VTtemp[fe+kk, jj] * O1 * (bbb * aa - a)
            VTtemp[fe+kk, jj] *= bbb
          end
        else
          b = sqrt(Z1' * Z1)
          tau = Z1 / b
          bb = rand(Chisq(maximum(size(Z1))))
          bbb = b/bb
          Ztemp[kk][Z00] = bb * tau
          vert = setminus(1:dim, fe + kk) # -->> !! setminus not defined !!
          for jj in 1:VC[i]
            VTtemp[fe+kk, jj] *= bbb
          end
        end
      end
      for ii in 1:ree
        ZZ[ii] = hcat(ZZ[ii], Ztemp[ii])
      end
      VTVT[i] = VTtemp
    end
    Z = ZZ
    VT = VTVT

    #----------------------------------------------------flip negatives
    for i in 1:N
      for j in 1:ree #only look at random effects
        if signs[j, i] == -1
          VT[i][fe+j, :] *= -1  #only need to flip the negatives
          Z[j][:, i] *= -1
        end
      end
    end
    #pick the coordinates
    VT_end = zeros(Float64, dim, N)
    for i in 1:N
      for j in 1:dim
        if rand() <= 0.5
          if j <= fe
            VT_end[j, i] = minimum(VT[i][j, :])
          else
            VT_end[j, i] = max(minimum(VT[i][j, :]), 0)
          end
        else
          if j <= fe
            VT_end[j, i] = maximum(VT[i][j, :])
          else
            VT_end[j, i] = max(maximum(VT[i][j, :]), 0)
          end
        end
      end
    end
    if k_n == K_n #if finished pick coordinates
      #pick the coordinates
      VT_end = zeros(Float64, dim, N)
      for i in 1:N
        for j = 1:dim
          if rand() <= 0.5
            if j <= fe
              VT_end[j, i] = minimum(VT[i][j, :])
            else
              VT_end[j, i] = max(minimum(VT[i][j, :]), 0)
            end
          else
            if j <= fe
              VT_end[j, i] = maximum(VT[i][j, :])
            else
              VT_end[j, i] = max(maximum(VT[i][j, :]), 0)
            end
          end
        end
      end
      VERTEX = VT_end
      WEIGHT = WT
      p = fe
      r = ree
    end
  end
  return (VERTEX, WEIGHT)
end





data = [2.0 3; 2 3; 3 4; 4 5; 4 5; 6 7]
FE = hcat([1.0; 1; 1; 1; 1; 1])
RE = hcat([1; 1; 2; 2; 3; 3])
N = 4

(VERTEX, WEIGHT) = gfilmm(data, FE, RE, 5000, 2500.0)

fname = "/home/stla/Work/R/gfilmm/inst/julia/weight.bin"; # random temporary filename
open(fname,"w") do f
  # Make sure we write 64bit integer in little-endian byte order
  write(f,WEIGHT)
end

fname = "/home/stla/Work/R/gfilmm/inst/julia/vertex.bin"; # random temporary filename
open(fname,"w") do f
  # Make sure we write 64bit integer in little-endian byte order
  write(f,VERTEX)
end
