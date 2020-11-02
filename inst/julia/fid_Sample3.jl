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
      if all(VTsum .>= U) || all(VTsum .<= L) #do not satisfy constaints
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

VT2 = -[1.0; 2.0; 3.0; 4.0]
VTsum = -[4.0; 3.0; 2.0; 1.0]
L = -4.5
U = -1.5
fid_sample(VT2, VTsum, U, L)
