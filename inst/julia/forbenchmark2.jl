using LinearAlgebra
using StaticArrays

function rchi(n)
   return norm(randn(n))
end

## sort array M with respect to colon j
sortbycol = function (M, j)
  local col, perm
  @inbounds col = M[:, j]
  perm = sortperm(col)
  return M[perm, :]
end

function sampleInt(n)
  return convert(Vector{Int64}, sortbycol(hcat(1:n, rand(n)), 2)[:,1])
end

function ff(cdd, VTsum, x, vtsumw, VTi, vtiw)
  lambda = (x - vtsumw) / (VTsum[cdd] - vtsumw)
  return lambda * VTi[:, cdd] + (1.0 - lambda) * vtiw
end

function fid_vertex(
  VTi::Array{R},
  CCi::Array{Int64},
  VTsum::Vector{R},
  Uk::R,
  Lk::R,
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
  VTtemp = zeros(R, dim, 0)
  vert = 0
  # check lower constraints first
  if lcheckl != 0
    @inbounds CA = CCi[:, checkl]
    @inbounds CB = CCi[:, whichl]
    INT = falses(2 * n, lcheckl)
    for ll in 1:lcheckl
      @inbounds INT[CA[:, ll], ll] .= 1
    end
    #
    @inbounds VTsum_wl = VTsum[whichl]
    @inbounds VTi_wl = VTi[:, whichl]
    #
    for ii in 1:lwhichl
      @inbounds INT2 = INT[CB[:, ii], :]
      @inbounds use = findall(count(==(1), INT2, dims=1)[1, :] .== dim-1)
      if length(use) > 0
        vert = vert + length(use)
        @inbounds inters = map(dd -> [CB[findall(INT2[:, dd] .== 1), ii]; k+n], use)
        CCtemp = hcat(CCtemp, hcat(inters...))
        @inbounds VTs = map(c -> ff(c, VTsum, Lk, VTsum_wl[ii], VTi, VTi_wl[:, ii]), checkl[use])
        VTtemp = hcat(VTtemp, hcat(VTs...))
      end
    end
  end
  # check upper constraints now
  if lchecku != 0
    @inbounds CA = CCi[:, checku]
    @inbounds CB = CCi[:, whichu]
    INT = falses(2 * n, lchecku)
    for ll in 1:lchecku
      @inbounds INT[CA[:, ll], ll] .= 1
    end
    #
    @inbounds VTsum_wu = VTsum[whichu]
    @inbounds VTi_wu = VTi[:, whichu]
    #
    for ii in 1:lwhichu
      @inbounds INT2 = INT[CB[:, ii], :]
      @inbounds use = findall(count(==(1), INT2, dims=1)[1, :] .== dim-1)
      if length(use) > 0
        vert = vert + length(use)
        @inbounds inters = map(dd -> [CB[findall(INT2[:, dd] .== 1), ii]; k], use)
        CCtemp = hcat(CCtemp, hcat(inters...))
        @inbounds VTs = map(c -> ff(c, VTsum, Uk, VTsum_wu[ii], VTi, VTi_wu[:, ii]), checku[use])
        VTtemp = hcat(VTtemp, hcat(VTs...))
      end
    end
  end
  if !isempty(both)
    @inbounds CCtemp = hcat(CCtemp, CCi[:, both])
    @inbounds VTtemp = hcat(VTtemp, VTi[:, both])
    vert = vert + length(both)
  end
  # bigs = VTtemp .> 5000
  # if any(bigs)
  #   println("bigs")
  # end
  return (VTtemp, CCtemp, vert)
end

function approx(x::R, n::Int64) where {R<:Real}
  return round(x * 10^n) / 10^n
end

function fid_sample(VT2::Vector{R}, VTsum::Vector{R}, U::R, L::R) where {R<:Real}
  high = findall(VT2 .> 0.0)
  low = findall(VT2 .< 0.0)
  zero = findall(VT2 .== 0.0)
  if (!isempty(low) && !isempty(high)) || !isempty(zero) # some positive and negative signs
    UU = map(sign, U .- VTsum)
    LL = map(sign, L .- VTsum)
    SS = map(sign, VT2)
    whichnot = findall(VT2 .!= 0)
    which = findall(VT2 .== 0)
    if length(zero) == length(VT2) # all are zero
      MAX = Inf
      MIN = -Inf
      temp = any(UU .== 1) * any(LL .== -1)
    elseif (
      !any(SS .== -1) && !any(UU[which] .== -1) && !any(LL[which] .== -1) # to improve c1 && d1 etc
    ) || (!any(SS .== 1) && !any(UU[which] .== 1) && !any(LL[which] .== 1))
      @inbounds Us = (U .- VTsum[whichnot]) ./ VT2[whichnot]
      @inbounds Ls = (L .- VTsum[whichnot]) ./ VT2[whichnot]
      MAX = Inf
      MIN = min(minimum(Us), minimum(Ls))
      temp = 0.5 - atan(MIN) / pi # weight adjustment
    elseif (
      !any(SS .== 1) && !any(UU[which] .== -1) && !any(LL[which] .== -1)
    ) || (
      !any(SS .== -1) && !any(UU[which] .== 1) && !any(LL[which] .== 1)
    )
      @inbounds Us = (U .- VTsum[whichnot]) ./ VT2[whichnot]
      @inbounds Ls = (L .- VTsum[whichnot]) ./ VT2[whichnot]
      MAX = max(maximum(Us), maximum(Ls))
      MIN = -Inf
      temp = atan(MAX) / pi + 0.5  # weight adjustment
    else
      if isempty(high)
        Hmax = -Inf
        Hmin = Inf
      else
        @inbounds HUs = (U .- VTsum[high]) ./ VT2[high]
        @inbounds HLs = (L .- VTsum[high]) ./ VT2[high]
        Hmax = max(maximum(HUs), maximum(HLs))
        Hmin = min(minimum(HUs), minimum(HLs))
      end
      if isempty(low)
        Lmax = -Inf
        Lmin = Inf
      else
        @inbounds LUs = (U .- VTsum[low]) ./ VT2[low]
        @inbounds LLs = (L .- VTsum[low]) ./ VT2[low]
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
        Pprob = 0.5 - atan(bottom_pos) / pi
      elseif bottom_pos == -Inf
        Pprob = atan(top_pos) / pi + 0.5
      end
      if top_neg == Inf
        Nprob = 0.5 - atan(bottom_neg) / pi
      elseif bottom_neg == -Inf
        Nprob = atan(top_neg) / pi + 0.5
      end
      temp = Pprob + Nprob
      Pprob /= temp
      Nprob = 1.0 - Pprob
      if rand() <= Nprob # choose negative range
        MIN = bottom_neg
        MAX = top_neg
      else # choose positive range
        MIN = bottom_pos
        MAX = top_pos
      end
    end
    #######################################################################
    y = atan(MAX) / pi
    x = atan(MIN) / pi
    u = x + 0.5 + (y - x) * rand()
    ZZ = tan(pi * (u - 0.5))  #Inverse cdf
    wt = exp(-ZZ^2 / 2.0) * (1.0 + ZZ^2) * temp
  else
    Us = (U .- VTsum) ./ VT2
    Ls = (L .- VTsum) ./ VT2
    MAX = max(maximum(Us), maximum(Ls))
    MIN = min(minimum(Us), minimum(Ls))
    y = atan(MAX) / pi
    x = atan(MIN) / pi
    u = x + 0.5 + (y - x) * rand()
    ZZ = tan(pi * (u - 0.5)) #Inverse cdf
    wt = exp(-ZZ^2 / 2.0) * (1.0 + ZZ^2) * (y - x)
  end
  return (ZZ, wt, MIN, MAX)
end


#srand(421)
#function gfilmm4(data::Array{R,2}, FE::Array{Float64,2}, RE::Array{Int64,2}, N::Int64, thresh::Float64) where {R<:Real}
function gfilmm4(data::Array{R,2}, FE, RE, N::Int64, thresh::Float64) where {R<:Real}

  local WT, VERTEX, WEIGHT, VTtemp, VCtemp

  n = size(data, 1)

  if size(data, 2) == 2
    L = data[:, 1]
    U = data[:, 2]
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

  RE = SMatrix{size(RE,1),size(RE,2)}(RE)
  FE = SMatrix{size(FE,1),size(FE,2)}(FE)

  dim = fe + ree  #dimension of the space
  rep = 1

  ########--------SET-UP ALGORITHM OBJECTS
  Z = Vector{Array{Float64,2}}(undef, ree)     #Particles
  weight = ones(Float64, n, N)
  ESS = N * ones(Float64, n)     #Effective sample size
  VC = zeros(Int64, N)       #Number of vertices
  VT = Vector{Array{Float64,2}}(undef, N)        #Verticies
  CC = Vector{Array{Int64,2}}(undef, N)        #Constraints

  ########--------SAMPLE ALL Z's

  A = Vector{Array{Float64,2}}(undef, N)
  for j in 1:N
    A[j] = @SMatrix zeros(Float64, n, 0)
  end

  for i in 1:ree  #Sample all the Z's
    Z[i] = randn(E[i], N)
    REblock = RE[:, (ESUM[i]-E[i]+1):ESUM[i]]
    # As = map(j -> REblock * Z[i][:, j], 1:N)
    # println(As)
    # A = hcat(As...)
    for j in 1:N
      A[j] = hcat(A[j], REblock * Z[i][:, j])
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
    @inbounds AA = hcat(vcat(FE, -FE), vcat(A[i], -A[i]))
    b = vcat(U, -L)
    for ii in 1:(2^dim) #number of vertices
      @inbounds II = USE[ii, :]  #constraints to use are by row
      @inbounds V[:, ii] = AA[II, :] \ b[II]
    end
    VT[i] = V
    CC[i] = USE'
  end

  VC = 2^dim * ones(Int64, N)

  ########--------MAIN ALGORITHM
  K_n = convert(Int64, ceil(length(K) / break_point))
  K_temp = Vector{Vector{Int64}}(undef, K_n)
  for i in 1:(K_n-1)
    @inbounds K_temp[i] = K[((i-1)*break_point+1):(i*break_point)]
  end
  @inbounds K_temp[K_n] = K[((K_n-1)*break_point+1):end]

  K1 = Vector{Int64}(undef, 0)
  for k_n in 1:K_n
    @inbounds K1 = vcat(K1, K_temp[k_n])
    for k in K_temp[k_n]
      level = k  #the level is the observation we are on
      if k_n > 1
        for i in 1:ree
          @inbounds Z[i] = vcat(Z[i], randn(E[i] - size(Z[i], 1), N))
        end
      end
      for i in 1:N
        @inbounds VT1 = VT[i]  #vertex values
        @inbounds VT2 = VT1[dim, :]  #value to be resampled
        # remove value to be resampled
        keepcols = setdiff(1:size(VT1, 1), dim)
        @inbounds VT1 = VT1[keepcols, :]  #
        if fe > 0
          @inbounds Z1 = FE[k, :]  #Assign fixed effect values
        else
          Z1 = zeros(Float64, 0)
        end
        for j in 1:ree
          @inbounds Z1 = vcat(Z1, Z[j][RE2[k, j], i])
        end
        @inbounds Z1 = Z1[keepcols]   #remove column of effect to be sampled
        VTsum = (Z1' * VT1)'
        @inbounds (Z[ree][k, i], weight[k, i]) = fid_sample(VT2, VTsum, U[k], L[k]) ###Sample -
        # if wt == 0.0
        #   println("wt 0")
        # end
        @inbounds VTsum += Z[ree][k, i] * VT2
        @inbounds VT1 = VT[i]
        @inbounds CC1 = CC[i]
        @inbounds (VT[i], CC[i], VC[i]) =
          fid_vertex(VT[i], CC[i], VTsum, U[k], L[k], dim, k, n)
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
          weight[k, i] = 0
        end
      end
      WT = vec(prod(weight, dims = 1)) #only last re is restricted
      if sum(WT) == 0#
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        println("!!!   Error: possible underflow. !!!")
        println("!!!   Potential resolutions:      !!!")
        println("!!!   1.  Widen data intervals   !!!")
        println("!!!   2.  Center data            !!!")
        println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        error("all weights are zero")
      end
      WT /= sum(WT)
      ESS[k] = 1 / (WT' * WT)
      #---------------------------------------------------------------Resample
      if ESS[k] < thresh && k < K[end]-1 #- 1 # -->> pourquoi -1 ?
        println("RESAMPLE")
        u = zeros(Float64, N)
        N_sons = zeros(Int64, N)
        # generate the cumulative distribution
        dist = N * cumsum(WT)
        aux = rand()    # sample uniform rv in [0 1]
        u = aux .+ (0:(N-1))
        j = 1
        for i in 1:N
          @inbounds while (u[i] > dist[j])
            j += 1
          end
          @inbounds N_sons[j] += 1
        end
        Nsons_sum = zeros(Int64, N+1)
        for i in 2:(N+1)
          @inbounds Nsons_sum[i] = Nsons_sum[i-1] + N_sons[i-1]
        end
        @inbounds Nsons_totalSum = Nsons_sum[N+1]
        JJ = union(1:k, K_start) # WRONG !!!
        ZZ = Vector{Array{Float64,2}}(undef, ree)
        for ii in 1:ree
          @inbounds ZZ[ii] = zeros(Float64, size(Z[ii], 1), Nsons_totalSum)
        end
        VCVC = zeros(Int64, N)
        CCCC = Vector{Array{Int64,2}}(undef, N)
        VTVT = Vector{Array{R,2}}(undef, N)
        for i in 1:N
          @inbounds if N_sons[i] > 0
            @inbounds VCtemp = VC[i] * ones(Int64, N_sons[i])
            Ztemp = Vector{Any}(undef, ree)
            VTtemp = Vector{Array{Float64,2}}(undef, N_sons[i])
            @inbounds copy = N_sons[i] - 1  # # to be resampled
            for ii in 1:N_sons[i]
              @inbounds VTtemp[ii] = VT[i]  #copy original vertices
            end
            for ii in 1:ree
              @inbounds Ztemp[ii] = SMatrix{size(Z[ii],1), N_sons[i]}(repeat(Z[ii][:, i], outer = (1, N_sons[i])))  #copy Z
            end
            ##########################################################
            if copy > 0
                ord = sampleInt(ree)
                for kk in ord
                  for ii in 1:copy
                    XX = @SMatrix zeros(Float64, n, 0)
                    for jj in 1:ree
                      @inbounds XX = hcat(XX, RE[:, (ESUM[jj]-E[jj]+1):ESUM[jj]] * Ztemp[jj][:, ii])
                    end
                    @inbounds XX = XX[JJ, setdiff(1:size(XX, 2), kk)] #remove column of effect to be resampled
                    @inbounds Z1 = Ztemp[kk][unique(RE2[JJ, kk]), ii]  #Z being resampled
                    @inbounds CO2 = RE[JJ, (ESUM[kk]-E[kk]+1):ESUM[kk]]
                    level0 = findall(vec(all(CO2 .== 0, dims=1)))
                    Z0 = findall(Z1 .== 0.0)  #Z's not sampled
                    Z00 = setdiff(1:length(Z1), Z0) #These are the levels with Z for effect kk
                    @inbounds CO2 = CO2[:, setdiff(1:size(CO2, 2), level0)]
                    @inbounds Z1 = Z1[Z00]#SVector{length(Z00)}(Z1[Z00])
                    if fe > 0
                      @inbounds XX = hcat(FE[JJ, :], XX)
                    end
                    MAT = hcat(-XX, CO2)
                    if rank(MAT) < size(MAT, 2)
                      NUL = nullspace(MAT) #
                      @inbounds n1 = SMatrix{size(NUL, 1) - size(CO2, 2), size(NUL, 2)}(NUL[1:size(NUL, 1) - size(CO2, 2), :])
                      @inbounds n2 = SMatrix{size(CO2, 2), size(NUL, 2)}(NUL[(size(NUL, 1) - size(CO2, 2)+1):end, :])
                      qrf = qr(n2)
                      @inbounds O2 = qrf.Q[:, 1:size(n2,2)]
                      RR = qrf.R
                      O1 = n1 * inv(UpperTriangular(RR))
                      a = O2' * Z1
                      O2a = O2 * a
                      M = Z1 - O2a
                      b = sqrt(M' * M)
                      tau = M / b
                      rankO2 = size(O2, 2)
                      bb = rchi(length(Z1) - rankO2)
                      bbb = b/bb #
                      aa = randn(rankO2, 1)
                      @inbounds Ztemp[kk][Z00, ii] = O2 * aa + bb * tau
                      vert = setdiff(1:dim, fe + kk)
                      @inbounds for jj in 1:VC[i]
                        @inbounds VTtemp[ii][vert, jj] -=
                          VTtemp[ii][fe+kk, jj] * O1 * (bbb * aa - a)
                        @inbounds VTtemp[ii][fe+kk, jj] *= bbb
                      end
                    else
                      b = sqrt(Z1' * Z1)
                      tau = Z1 / b
                      bb = rchi(length(Z1))
                      @inbounds Ztemp[kk][Z00, ii] = bb * tau
                      vert = setdiff(1:dim, fe + kk)
                      @inbounds for jj in 1:VC[i]
                        @inbounds VTtemp[ii][fe+kk, jj] *= b / bb
                      end
                    end
                    # bigs = VTtemp[ii] .> 5000
                    # if any(bigs)
                    #   println("bbigs")
                    # end
                  end
                end
            end
            for ii in 1:ree
              @inbounds ZZ[ii][:, (Nsons_sum[i]+1):Nsons_sum[i+1]] = Ztemp[ii]
            end
            @inbounds VCVC[(sum(N_sons[1:(i-1)])+1):sum(N_sons[1:i])] = VCtemp
            @inbounds d = Nsons_sum[i+1]#sum(N_sons[1:(i-1)]) # i in C++! why i-1 here ?
            @inbounds for kk in 1:N_sons[i]
              @inbounds VTVT[kk+d] = VTtemp[kk]
              @inbounds CCCC[kk+d] = CC[i]
            end
          end
        end
        Z = ZZ
        VT = VTVT
        VC = VCVC
        CC = CCCC
        weight = ones(Float64, n, N)  #assign weights of error matrix to 1
      end
    end  #ends resampling for k=K1

    #----------------------------------------------------determine signs
    signs = zeros(Int64, ree, N)
    for i in 1:N
      for j in 1:ree
        @inbounds positive = VT[i][fe+j, :] .> 0
        @inbounds negative = VT[i][fe+j, :] .< 0
        if count(positive) == VC[i] #i.e. all are positive
          @inbounds signs[j, i] = 1
        elseif count(negative) == VC[i] #i.e. all are negative
          @inbounds signs[j, i] = -1
        end
      end
    end

    #----------------------------------------------------FINAL RESAMPLE
    ZZ = Vector{Array{Float64,2}}(undef, ree)
    VTVT = Vector{Array{Float64,2}}(undef, N)

    n1 = sort(vcat(K1, K_start))
    nn = Vector{Vector{Int64}}(undef, ree)
    for ii in 1:ree
      @inbounds nn[ii] = unique(RE2[n1, ii])
      @inbounds ZZ[ii] = @SMatrix zeros(Float64, length(nn[ii]), N)
    end

    for i in 1:N
      Ztemp = Vector{Any}(undef, ree)
      @inbounds VTtemp = VT[i]
      for ii in 1:ree
        @inbounds Ztemp[ii] = Z[ii][nn[ii], i]  #copy Z
      end
      ord = sampleInt(ree)
      for kk in ord
        XX = zeros(Float64, size(RE, 1), 0)
        eff = setdiff(1:ree, kk)
        for jj in eff
          @inbounds XX = hcat(
            XX, RE[:, (ESUM[jj]-E[jj]+1):(ESUM[jj]-E[jj]+length(nn[jj]))]
          ) * Ztemp[jj]
        end
        @inbounds Z1 = Ztemp[kk]  #Z being resampled
        @inbounds CO2 = RE[:, (ESUM[kk]-E[kk]+1):(ESUM[kk]-E[kk]+length(nn[kk]))]
        Z0 = findall(Z1 .== 0)  #Z's not sampled
        Z00 = setdiff(1:length(Z1), Z0) #These are the levels with Z for effect kk
        @inbounds Z1 = SVector{length(Z00)}(Z1[Z00])
        if fe > 0
          XX = hcat(FE, XX)
        end
        MAT = hcat(-XX, CO2)
        if rank(MAT) < size(MAT, 2)
          NUL = nullspace(MAT)
          @inbounds n1 = SMatrix{size(NUL, 1) - size(CO2, 2), size(NUL, 2)}(NUL[1:size(NUL, 1) - size(CO2, 2), :])
          @inbounds n2 = SMatrix{size(CO2, 2), size(NUL, 2)}(NUL[(size(NUL, 1) - size(CO2, 2)+1):end, :])
          qrf = qr(n2)
          @inbounds O2 = qrf.Q[:, 1:size(n2, 2)]
          RR = qrf.R
          O1 = n1 * inv(UpperTriangular(RR))
          a = O2' * Z1
          O2a = O2 * a
          M = Z1 - O2a
          b = sqrt(M' * M)
          tau = M / b
          rankO2 = size(O2, 2)
          bb = rchi(length(Z1) - rankO2)
          bbb = b/bb
          aa = randn(rankO2, 1)
          @inbounds Ztemp[kk][Z00] = O2 * aa + bb * tau
          vert = setdiff(1:dim, fe + kk)
          @inbounds for jj = 1:VC[i]
            @inbounds VTtemp[vert, jj] -= VTtemp[fe+kk, jj] * O1 * (bbb * aa - a)
            @inbounds VTtemp[fe+kk, jj] *= bbb
          end
        else
          b = sqrt(Z1' * Z1)
          tau = Z1 / b
          bb = rchi(length(Z1))
          bbb = b/bb
          @inbounds Ztemp[kk][Z00] = bb * tau
          @inbounds for jj in 1:VC[i]
            @inbounds VTtemp[fe+kk, jj] *= bbb
          end
        end
        # bigs = VTtemp .> 5000
        # if any(bigs)
        #   println("bbbigs")
        # end
      end
      for ii in 1:ree
        @inbounds ZZ[ii][:,i] = Ztemp[ii]
      end
      @inbounds VTVT[i] = VTtemp
    end
    Z = ZZ
    VT = VTVT

    #----------------------------------------------------flip negatives
    for i in 1:N
      for j in 1:ree #only look at random effects
        @inbounds if signs[j, i] == -1
          @inbounds VT[i][fe+j, :] *= -1  #only need to flip the negatives
          @inbounds Z[j][:, i] *= -1
        end
      end
    end

    if k_n == K_n #if finished pick coordinates
      #pick the coordinates
      VT_end = zeros(Float64, dim, N)
      for i in 1:N
        for j = 1:dim
          if rand() < 0.5
            if j <= fe
              @inbounds VT_end[j, i] = minimum(VT[i][j, :])
            else
              @inbounds VT_end[j, i] = max(minimum(VT[i][j, :]), 0)
            end
          else
            if j <= fe
              @inbounds VT_end[j, i] = maximum(VT[i][j, :])
            else
              @inbounds VT_end[j, i] = max(maximum(VT[i][j, :]), 0)
            end
          end
        end
      end
      VERTEX = VT_end
      WEIGHT = WT
    end
  end
  return (VERTEX, WEIGHT, ESS)
end

data = [2.0 3; 2 3; 3 4; 4 5; 4 5; 6 7]
FE = hcat([1.0; 1; 1; 1; 1; 1])
RE = hcat([1; 1; 2; 2; 3; 3])
N = 4

function f(N)
  (VERTEX, WEIGHT, ESS) = gfilmm4(data, FE, RE, N, convert(Float64, N/2))
  return VERTEX
end


# function test1()
#   M0 = ones(10, 10)
#   M = ones(10,0)
#   for i in 1:100
#     M = hcat(M,M0)
#   end
#   return M
# end
#
# function test2()
#   M = ones(10, 10)
#   Ms = map(i -> M, 1:100)
#   out = hcat(Ms...)
#   return out
# end
#
# using BenchmarkTools
#
# @btime test1()
# @btime test2()

using LinearAlgebra

function f1(M)
  x = qr!(M)
  Q = x.Q;
  R = x.R;
  return size(Q,1) + size(R,1)
end

function f2(M)
  Q, R = qr!(M)
  return size(Q,1) + size(R,1)
end

using BenchmarkTools

M = ones(10, 10)

@btime f1(M)

@btime f2(M)
