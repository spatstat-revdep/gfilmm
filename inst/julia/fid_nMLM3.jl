using LinearAlgebra
using Distributions # to sample chi²

AAA = []
AAAA = []
KKK = []
OOO2 = []
ZZZ1 = []
COOO2 = []
MMMAT = []
ZZZZZ = []
iiii = []
ZZZZ = cell(16)
ESSS = cell(16)
srand(421)

data = [2.0 3; 3 4; 3 4; 4 5; 5 6; 6 7]
FE = [1.0; 1; 1; 1; 1; 1]
RE = [1; 1; 1; 2; 2; 2]
N = 5

fid_nMLM = function (data::Array{R,2}, FE::Array{R,2}, RE::Array{Int64,2}, N::Int64, thresh::R) where {R<:Real}
  global AAA, AAAA, KKK, OOO2, COOO2, ZZZ1, MMMAT, ESSS, ZZZZZ, iiii, ZZZZ
  VERTEX = []
  WEIGHT = []
  ESS = []
  p = []
  r = []

  n = size(data,1)
  random_design = 1
  Y = data
  if size(data,2) == 2
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
      A_temp = vcat(AT, AA[k, :])
      if rank(A_temp) > r
        AT = A_temp
        r += 1
        C = vcat(C, k)
      end
    end
  end

  K = setdiff(1:n, C)
  K_start = C
  Z[ree][K, :] = 0  #remove error particles not selected in intitialization
  # K will be the values left to be sampled

  ########--------FIND INITIAL VERTICES

  USE = repeat(C, outer = (2^dim, 1))  #all combinations of I
  for j in 1:dim
    if j == 1
      USE[(2^(dim-1)+1):(2^dim), 1] = USE[1:(2^(dim-1)), 1] + n
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
      b = vcat(U, -L)
      AA = AA[II, :]
      b = b[II]
      V[:, ii] = AA \ b
    end
    VT[i] = V
  end
  for k in 1:N
    CC[k] = reshape(USE', dim, 2^dim)  #constraints are the same for all N particles
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
    K1 = hcat(K1, K_temp[k_n])
    for k in K_temp[k_n]
      effect = ree  #only have the error term left
      level = k  #the level is the observation we are on
      if k_n > 1
        for i in 1:ree
          Z[i] = vcat(Z[i], randn(E[i] - size(Z[i], 1), N))
        end
      end
      for i in 1:N
#        m = convert(Int, VC[1, i])
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
        for j = 1:ree
          if RE2[k, j] == 0 # -->> it's never 0, I think
            Z1 = [Z1 0]
          else
            Z1 = hcat(Z1, Z[j][RE2[k, j], i])
          end
        end
        Z1 = Z1[keepcols]   #remove column of effect to be sampled
        VTsum = (Z1' * VT1)'
        (ZZ, wt) = fid_sample(VT2, VTsum, U[k], L[k]) ###Sample -
        Z[effect][k, i] = ZZ
        weight[effect][k, i] = wt # -->> only one 'weight' object is used !
        VTsum += Z[effect][k, i] * VT2
        VT1 = VT[i]
        CC1 = CC[i]
#        (VTtemp, CCtemp, vert) = fid_vertex(VT1, CC1, VTsum, U[k], L[k], dim, k, n) # -->> I removed the 'm'
#        VC[i] = vert
#        CC[i] = CCtemp
#        VT[i] = VTtemp
        (VT[i], CC[i], VC[i]) = fid_vertex(VT1, CC1, VTsum, U[k], L[k], dim, k, n) # -->> I removed the 'm'
        if VC[i] == 0 #check
          weight[effect][k, i] = 0
        end
      end
      # ZZZZ[k] = Z[effect, 1] -->> not used
      WT = cumprod(weight[ree]) #only last re is restricted
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
      ESS[k] = 1 / norm(WT)^2
      # ESSS[k] = ESS[k, 1] -->> notused
      #---------------------------------------------------------------Resample
      if ESS[k] < thresh && k < K[end] - 1
        u = zeros(N)
        N_sons = zeros(Int, N)
        # generate the cumulative distribution
        dist = N * cumsum(WT)
        aux = rand(1)    # sample uniform rv in [0 1]
        u = aux + (0:(N-1))
        #u = u ./ N -->> I multiplied 'dist' by N instead
        j = 1
        for i in 1:N
          while (u[i] > dist[j])
            j += 1
          end
          N_sons[j] += 1
        end
        #KEEP[k,:]=N_sons
        II = setdiff(1:n, union(1:k, C))
        JJ = setdiff(1:n, II)
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
                ord = sortbycol([[1:ree]'; rand(1, ree)]', 2)
                ord = ord[:, 1]'  #Order to resample.  Each re will be resampled.
                ord = convert(Array{Int,1}, vec(ord))
                for kk in ord
                  for ii in 1:copy
                    CO = RE[JJ, :]
                    XX = zeros(Float64, n, 0)
                    for jj in 1:ree
                      XX = hcat(XX, RE[:, (ESUM[jj]-E[jj]+1):ESUM[jj]] * Ztemp[jj][:, ii])
                    end
                    XX = XX[JJ, :]
                    XX = XX[:, setdiff(1:size(XX, 2), kk)] #remove column of effect to be resampled
                    temp = RE2[JJ, kk] .!= 0  #find which levels of kk have been sampled
                    Z1 = Ztemp[kk][unique(RE2[JJ[temp], kk]), ii]  #Z being resampled
                    CO2 = RE[JJ, [(ESUM[kk]-E[kk]+1):ESUM[kk]]]
                    level0 = findall(sum(abs(CO2), 1) .== 0) #levels not sampled yet # -->> TO REPLACE WITH count
                    Z0 = findall(Z1 .== 0.0)  #Z's not sampled
                    Z00 = setdiff(1:length(Z1), Z0) #These are the levels with Z for effect kk
                    CO2 = CO2[:, setminus(1:size(CO2, 2), level0)]
                    Z1 = setdiff(Z1, Z0) # -->> was setminus(Z0, Z1) --- maybe use setdiff!
                    if fe > 0
                      XX = hcat(FE[JJ, :], XX)
                    end
                    MAT = hcat(-XX, CO2)
                    if rank(MAT) < size(MAT, 2)
                      NUL = nullspace(MAT) #
                      n1 = NUL[1:size(NUL, 1) - size(CO2, 2), :]
                      n2 = NUL[(size(NUL, 1) - size(CO2, 2)+1):end, :]
                      O2 = qr(n2).Q
                      B = CO2 * O2
                      O1 = XX \ B
                      a = O2' * Z1
                      b = sqrt((Z1 - O2 * a)' * (Z1 - O2 * a))
                      tau = (Z1 - O2 * a) / b
                      AAA = maximum(size(Z1)) - rank(O2)
                      bb = rand(Chisq(maximum(size(Z1)) - rank(O2))) # -->> rank O2 is ncol or nrow
#                      bb = sqrt(randchi2(max(size(Z1)) - rank(O2)))
                      bbb = (b/bb)[1] #
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
                      bb = rand(Chisq(maximum(size(Z1))))
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
              ZZ[ii, 1] = hcat(ZZ[ii, 1], Ztemp[ii, 1])  # ? c'est pareil que ZZ[ii,1]=2temp[ii,1]
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
    ZZ = cell(ree, 1)
    VTVT = cell(1, N)
    #ZZZZZ = Z -->> not used
    for i = 1:N
      iiii = i
      Ztemp = cell(ree, 1)
      VTtemp = VT[1, i]
      n1 = sort([K1; K_start])
      nn = cell(ree, 1)
      for ii = 1:ree
        nn[ii, 1] = unique(RE2[vec(n1), ii])
        Ztemp[ii, 1] = Z[ii, 1][nn[ii, 1], i]  #copy Z
        if i == 1 # initialisation que pour i=1
          ZZ[ii, 1] = zeros(size(Ztemp[ii, 1], 1), 0)
        end
      end
      ord = sortbycol([[1:ree]'; rand(1, ree)]', 2)
      ord = ord[:, 1]'  #Order to resample.  Each re will be resampled.
      ord = convert(Array{Int,2}, ord)
      for kk in ord
        CO = RE
        XX = zeros(size(RE, 1), 0)
        eff = [1:ree]
        eff = setminus(kk, eff)
        for jj in eff
          XX = [XX RE[:, ESUM[jj]-E[jj]+1:ESUM[jj]-E[jj]+length(nn[jj, 1])] *
              Ztemp[jj, 1]]
        end
        #temp=find(RE2[:,kk].!=0)  #find which levels of kk have been sampled
        Z1 = Ztemp[kk, 1]  #Z being resampled
        CO2 = RE[:, ESUM[kk]-E[kk]+1:ESUM[kk]-E[kk]+length(nn[kk, 1])]  # c'?tait E[:,kk] mais une seule ligne
        level0 = find(sum(abs(CO2), 1) .== 0) #levels not sampled yet
        Z0 = find(Z1 .== 0)  #Z's not sampled
        Z00 = [1:length(Z1)]
        Z00 = setminus(Z0, Z00) #These are the levels with Z for effect kk
        CO2 = CO2[:, setminus(level0, 1:size(CO2, 2))]
        Z1 = setminus(Z0, Z1)
        if fe > 0
          XX = [FE XX]
        end
        MAT = [-XX CO2]
        MMMAT = MAT
        COOO2 = CO2
        if rank(MAT) < length(MAT[1, :])
          NUL = null(MAT) #
          n1 = NUL[1:length(NUL[:, 1])-length(CO2[1, :]), :]
          n2 = NUL[length(NUL[:, 1])-length(CO2[1, :])+1:end, :]
          O2 = orth(n2)  # rank(O2) = ncol(NUL) ?
          B = CO2 * O2
          O1 = XX \ B
          a = O2' * Z1
          b = sqrt((Z1 - O2 * a)' * (Z1 - O2 * a))
          tau = (Z1 - O2 * a) / b
          OOO2 = O2
          ZZZ1 = Z1
          AAAA = max(size(Z1)) - rank(O2)
          bb = sqrt(randchi2(max(size(Z1)) - rank(O2)))
          bbb = (b/bb)[1] #
          aa = randn(rank(O2), 1)
          Ztemp[kk, 1][Z00, 1] = O2 * aa + bb * tau
          vert = setminus(fe + kk, [1:dim])
          for jj = 1:VC[i]
            VTtemp[vert, jj] = VTtemp[vert, jj] - VTtemp[fe+kk, jj] * O1 * (bbb * aa - a)
            VTtemp[fe+kk, jj] = (VTtemp[fe+kk, jj]*b/bb)[1]
          end
        else
          b = sqrt((Z1)' * (Z1))
          tau = (Z1) / b
          bb = sqrt(randchi2(max(size(Z1))))
          Ztemp[kk, 1][Z00, 1] = bb * tau
          vert = setminus(fe + kk, [1:dim])
          for jj = 1:VC[i]
            VTtemp[fe+kk, jj] = (VTtemp[fe+kk, jj]*b/bb)[1]
          end
        end
      end
      for ii = 1:ree
        ZZ[ii, 1] = [ZZ[ii, 1] Ztemp[ii, 1]]
      end
      VTVT[1, i] = VTtemp
    end
    Z = ZZ
    VT = VTVT
    #----------------------------------------------------flip negatives
    for i = 1:N
      for j = 1:ree #only look at random effects
        if signs[j, i] == -1
          VT[1, i][fe+j, :] = -VT[1, i][fe+j, :]  #only need to flip the negatives
          Z[j, 1][:, i] = -Z[j, 1][:, i]
        end
      end
    end
    #pick the coordinates
    VT_end = zeros(fe + ree, N)
    for i = 1:N
      for j = 1:ree+fe
        if rand(1)[1] <= 0.5
          if j <= fe
            VT_end[j, i] = min(VT[1, i][j, :])
          else
            VT_end[j, i] = max(min(VT[1, i][j, :]), 0)
          end
        else
          if j <= fe
            VT_end[j, i] = max(VT[1, i][j, :])
          else
            VT_end[j, i] = max(max(VT[1, i][j, :]), 0)
          end
        end
      end
    end
    if k_n == K_n #if finished pick coordinates
      #pick the coordinates
      VT_end = zeros(fe + ree, N)
      for i = 1:N
        for j = 1:ree+fe
          if rand(1)[1] <= 0.5
            if j <= fe
              VT_end[j, i] = min(VT[1, i][j, :])
            else
              VT_end[j, i] = max(min(VT[1, i][j, :]), 0)
            end
          else
            if j <= fe
              VT_end[j, i] = max(VT[1, i][j, :])
            else
              VT_end[j, i] = max(max(VT[1, i][j, :]), 0)
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
