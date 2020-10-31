function fid_sample(VT2, VTsum, U, L)
  high = find(VT2 .> 0)
  low = find(VT2 .< 0)
  zero = find(VT2 .== 0)
  MIN = []
  MAX = []
  bottom_pos = []
  top_pos = []
  bottom_neg = []
  top_neg = []
  ZZ = []
  wt = []
  temp = []
  if (isempty(low) == 0 && isempty(high) == 0) || (isempty(zero) == 0) #some positive and negative sigs
    UU = sign(U - VTsum)
    LL = sign(L - VTsum)
    SS = sign(VT2)
    whichnot = find(VT2 .!= 0)
    which = find(VT2 .== 0)
    if length(zero) == length(VT2) #all are zero
      if (all(VTsum .>= U)) || (all(VTsum .<= L)) #do not satisfy constaints
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
      !any(SS .== -1) == 1 && !any(UU[which] .== -1) == 1 && !any(LL[which] .== -1) == 1
    ) || (!any(SS .== 1) == 1 && !any(UU[which] .== 1) == 1 && !any(LL[which] .== 1) == 1)
      Us = (U - VTsum[whichnot]) ./ VT2[whichnot]'
      Ls = (L - VTsum[whichnot]) ./ VT2[whichnot]'
      MAX = Inf
      MIN = min(min(Us), min(Ls))
      temp = 1 - (atan(MIN) / pi + 0.5)  #weight adjustment
    elseif (
      !any(SS .== 1) == 1 && !any(UU[which] .== -1) == 1 && !any(LL[which] .== -1) == 1
    ) || (
      !any(SS .== -1) == 1 && !any(UU[which] .== 1) == 1 && !any(LL[which] .== 1) == 1
    )
      Us = (U - VTsum[whichnot]) ./ VT2[whichnot]'
      Ls = (L - VTsum[whichnot]) ./ VT2[whichnot]'
      MAX = max(max(Us), max(Ls))
      MIN = -Inf
      temp = atan(MAX) / pi + 0.5  #weight adjustment
    else
      HUs = (U - VTsum[high]) ./ VT2[high]'
      HLs = (L - VTsum[high]) ./ VT2[high]'
      Hmax = max(max(HUs), max(HLs))
      Hmin = min(min(HUs), min(HLs))
      LUs = (U - VTsum[low]) ./ VT2[low]'
      LLs = (L - VTsum[low]) ./ VT2[low]'
      Lmax = max(max(LUs), max(LLs))
      Lmin = min(min(LUs), min(LLs))
      if 0 <= roundn(Lmin - Hmax, -12)    ## roundn indisponible 
        bottom_pos = -Inf
        top_pos = Hmax
        bottom_neg = Lmin
        top_neg = Inf
      elseif roundn(Hmin - Lmax, -12) >= 0
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
        Pprob = 1 - (atan(bottom_pos) / pi + 0.5)
      elseif bottom_pos == -Inf
        Pprob = atan(top_pos) / pi + 0.5
      end
      if top_neg == Inf
        Nprob = 1 - (atan(bottom_neg) / pi + 0.5)
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
    wt = exp(-ZZ^2 / 2) * (1 + ZZ^2) * temp
  else
    Us = (U - VTsum) ./ VT2'
    Ls = (L - VTsum) ./ VT2'
    MAX = max(max(Us), max(Ls))
    MIN = min(min(Us), min(Ls))
    y = atan(MAX) / pi + 0.5 #cdf
    x = atan(MIN) / pi + 0.5
    u = x + (y - x) * rand()
    ZZ = tan(pi * (u - 0.5)) #Inverse cdf
    wt = exp(-ZZ^2 / 2) * (1 + ZZ^2) * (y - x)
  end
  return (ZZ, wt, MIN, MAX)
end