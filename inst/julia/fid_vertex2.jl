function fid_vertex(VT1, CC1, VTsum, U, L, m, dim, k, n)
  # VT1 : matric 4x16
  # CC1: idem integer entries
  # VTsum: colonne 16 
  # U,L : nbm
  # m,dim,k,n : integers
  whichVT = zeros(Int, m, 2)
  whichVT[:, 1] = VTsum .>= L
  whichVT[:, 2] = VTsum .<= U
  whichl = find(whichVT[:, 1] .== 1)  #vertices that satisfy lower constraint
  whichu = find(whichVT[:, 2] .== 1)  #vertices that satisfy upper constraint
  both = find(sum(whichVT, 2) .== 2) ## -->> findall(all(whichVT, dims=2))
                                     ##       -> cartesian index
                                     ##      findall(all(whichVT, dims=2)[:,1])
  checkl = find(whichVT[:, 1] .== 0) 
  checku = find(whichVT[:, 2] .== 0)
  CCtemp = zeros(dim, 0) # ou utiliser peut-?tre cell()
  VTtemp = zeros(dim, 0)
  vert = 0
  CA = CC1[:, checkl]
  CB = CC1[:, whichl]
  if isempty(checkl) == 0 #i.e. they do not all satisfy the lower constraint
    INT = zeros(Int, 2 * n, length(checkl))
    for ll = 1:length(checkl) #check lower constraints first
      INT[CA[:, ll], ll] = 1
    end
    for ii = 1:length(whichl)
      use = [] # local use  # ?
      INT2 = INT[CB[:, ii], :]
      use = find(sum(INT2, 1) .== dim - 1) # sum sur colonne
      for dd = 1:length(use)
        inter = CB[find(INT2[:, use[dd]] .== 1), ii]  #this will be intersection
        vert = vert + 1
        CCtemp = [CCtemp [inter; k + n]]  # ATTENTION A VERIFIER need to add n indicating lower constraint
        lambda = (L - VTsum[whichl[ii]]) / (VTsum[checkl[use[dd]]] - VTsum[whichl[ii]])
        VTtemp =
          [VTtemp lambda * VT1[:, checkl[use[dd]]] + (1 - lambda) * VT1[:, whichl[ii]]] # il faut coller en colonnes
      end
    end
  end
  CA = CC1[:, checku]
  CB = CC1[:, whichu]
  if isempty(checku) == 0 #i.e. they do not all satisfy the lower constraint
    INT = zeros(2 * n, length(checku))
    for ll = 1:length(checku) #check lower constraints first
      INT[CA[:, ll], ll] = 1
    end
    for ii = 1:length(whichu)
      use = []
      INT2 = INT[CB[:, ii], :]
      use = find(sum(INT2, 1) .== dim - 1)
      for dd = 1:length(use)
        inter = CB[find(INT2[:, use[dd]] .== 1), ii]   #this will be intersection
        vert = vert + 1
        CCtemp = [CCtemp [inter; k]]  #need to add n indicating lower constraint
        lambda = (U - VTsum[whichu[ii]]) / (VTsum[checku[use[dd]]] - VTsum[whichu[ii]])
        VTtemp =
          [VTtemp lambda * VT1[:, checku[use[dd]]] + (1 - lambda) * VT1[:, whichu[ii]]]
      end
    end
  end
  if isempty(both) == 0
    for ll in both'
      vert = vert + 1
      CCtemp = [CCtemp CC1[:, ll]]
      VTtemp = [VTtemp VT1[:, ll]]
    end
  end
  return (VTtemp, CCtemp, vert)
end