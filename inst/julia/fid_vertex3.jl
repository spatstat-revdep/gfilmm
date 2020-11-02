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
      print(INT2)
      use = findall(vec(sum(INT2, dims=1) .== dim - 1)) # sum sur colonne
      print(use)
      for dd = 1:length(use)
        inter = CB[findall(INT2[:, use[dd]] .== 1), ii]  #this will be intersection
        print(inter)
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
      print("\nINT2 below\n")
      INT2 = INT[CB[:, ii], :]
      use = findall((sum(INT2, dims=1) .== dim - 1)[1, :])
      print("\n************use*********:\n")
      print(use)
      print("\n")
      print(CB)
      for dd in 1:length(use)
        #print("\naaaaaaa\n")
        #print(use[dd])
        print("\n")
        #print(findall((INT2[:, use[dd]] .== 1)[:,1]))
        #print("\nxxx\n")
        inter = CB[findall((INT2[:, use[dd]] .== 1)[:,1]), ii]
        #print(inter)  #this will be intersection
        print("\n")
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
    print("\n!isempyboth\n")
    print(both)
    for ll in both # both' ??
      vert = vert + 1
      CCtemp = [CCtemp CC1[:, ll]]
      VTtemp = [VTtemp VT1[:, ll]]
    end
  end
  return (VTtemp, CCtemp, vert)
end

function fid_vertex2(
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
  m = length(VTsum) # ? pourquoi m en argument alors ?
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
  CCtemp = zeros(Int64, dim, 0) # ou utiliser peut-?tre cell()
  VTtemp = zeros(R, dim, 0)
  vert = 0
  CA = CC1[:, checkl]
  CB = CC1[:, whichl]
  if !isempty(checkl) #i.e. they do not all satisfy the lower constraint
    INT = zeros(Int64, 2 * n, length(checkl))
    INT[CA[:, 1:length(checkl)], 1:length(checkl)] .= 1
    for ii = 1:length(whichl)
      INT2 = INT[CB[:, ii], :]
      use = findall(sum(INT2, dims=1) .== dim - 1) # sum sur colonne
      for dd = 1:length(use)
        inter = CB[findall(INT2[:, use[dd]] .== 1), ii]  #this will be intersection
        vert = vert + 1
        CCtemp = hcat(CCtemp, [inter; k + n])  # ATTENTION A VERIFIER need to add n indicating lower constraint
        lambda = (L - VTsum[whichl[ii]]) /
          (VTsum[checkl[use[dd]]] - VTsum[whichl[ii]])
        VTtemp = hcat(
          VTtemp,
          lambda * VT1[:, checkl[use[dd]]] + (1.0 - lambda) * VT1[:, whichl[ii]]
        )
      end
    end
  end
  CA = CC1[:, checku]
  CB = CC1[:, whichu]
  if !isempty(checku) #i.e. they do not all satisfy the lower constraint
    INT = zeros(2 * n, length(checku))
    INT[CA[:, 1:length(checku)], 1:length(checku)] .= 1
#    for ll = 1:length(checku) #check lower constraints first
#      INT[CA[:, ll], ll] .= 1
#    end
    for ii = 1:length(whichu)
      INT2 = INT[CB[:, ii], :]
      use = findall(sum(INT2, dims=1) .== dim - 1)
      for dd = 1:length(use)
        inter = CB[findall(INT2[:, use[dd]] .== 1), ii]   #this will be intersection
        vert = vert + 1
        CCtemp = hcat(CCtemp, [inter; k])  #need to add n indicating lower constraint
        lambda = (U - VTsum[whichu[ii]]) /
          (VTsum[checku[use[dd]]] - VTsum[whichu[ii]])
        VTtemp = hcat(
          VTtemp,
          lambda * VT1[:, checku[use[dd]]] + (1.0 - lambda) * VT1[:, whichu[ii]]
        )
      end
    end
  end
  if !isempty(both)
    CCtemp = hcat(CCtemp, CC1[:, both])
    VTtemp = hcat(VTtemp, VT1[:, both])
    vert = vert + length(both)
  end
  return (VTtemp, CCtemp, vert)
end


VT1 = [
  1.0 2 3 4 5 6 7 8;
  8 7 6 5 4 3 2 1;
  5 5 6 7 7 4 5 4
]
CC1 = [
  3 3 2 2 1 1 2 2; # plante si head = 3 3
  2 1 2 1 2 1 2 1;
  1 1 1 1 2 2 2 2
]
VTsum = [1.0; 2; 4; 5; 9; 3; 5; 7]
U = 6.0
L = 2.0
m = 8
dim = 3
k = 2
n = 18

fid_vertex(VT1, CC1, VTsum, U, L, m, dim, k, n)
fid_vertex2(VT1, CC1, VTsum, U, L, m, dim, k, n)


2.0  3.0  4.0  6.0  7.0
7.0  6.0  5.0  3.0  2.0
5.0  6.0  7.0  4.0  5.0
,
1  2  2  1  2
1  2  1  1  2
1  1  1  2  2
,
5

#
3.5   3.71429  2.0  3.0  4.0  6.0  7.0
5.5   5.28571  7.0  6.0  5.0  3.0  2.0
6.25  6.14286  5.0  6.0  7.0  4.0  5.0
,
2.0  1.0  3.0  2.0  2.0  1.0  2.0
 1.0  1.0  1.0  2.0  1.0  1.0  2.0
 2.0  2.0  1.0  1.0  1.0  2.0  2.0
,
7
