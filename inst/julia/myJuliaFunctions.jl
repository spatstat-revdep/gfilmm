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
## remove the elements of B belonging to A e.g. setminus0([2 3 4], [1:7]) = [1 5 6 7]
# -->> setdiff([1:7], [2 3 4])
function setminus0(A, B)
  BminusA = []
  for k = 1:length(B)
    if !contains(A, B[k])
      BminusA = [BminusA, B[k]]
    end
  end
  return (BminusA)
end
## remove the elements of B whose indices belong to A
## e.g. setminus([2 3 4], ["one" "two" "three" "four" "five" "six" "seven"]) = ["one" "five" "six" "seven"]
function setminus(A, B)
  return (B[setminus0(A, [1:length(B)])])
end
## similiar to R or Matlab unique() function
## -->> it is in Julia!
function unique(x)
  local n
  n = length(x)
  y = []
  for i = 1:n
    local c
    c = count(x[i:n] .== x[i])
    if c == 1
      y = [y, x[i]]
    end
  end
  return (y)
end
