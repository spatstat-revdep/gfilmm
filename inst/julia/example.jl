#load("myJuliaFunctions.jl")
#load("fid_Sample2.jl")
#load("fid_vertex2.jl")
#load("fid_nMLM2.jl")
data = [
  1.624 1.625
  2.209 2.210
  2.097 2.098
  0.558 0.559
  -0.335 -0.334
  -0.971 -0.970
  -1.650 -1.649
  -2.338 -2.337
  -3.290 -3.289
  -4.291 -4.290
  2.862 2.863
  2.023 2.024
  -2.336 -2.335
  -0.613 -0.612
  -0.907 -0.906
  0.354 0.355
]
FE = ones(length(data[:, 1]), 1)
RE = [1 1; 1 1; 1 1; 1 1; 1 2; 1 2; 1 2; 1 2; 2 3; 2 3; 3 4; 3 4; 4 5; 4 5; 5 6; 5 6]
N = 1000
thresh = 500

(VERTEX, WEIGHT) = fid_nMLM(data, FE, RE, N, thresh)
#print(VERTEX)
#print(WEIGHT)


#writecsv("juliavertex.csv",VERTEX)
#writecsv("juliaweight.csv",WEIGHT)