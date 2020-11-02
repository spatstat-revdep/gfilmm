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
      println(use)
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
