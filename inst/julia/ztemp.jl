function f()
	local x # crash without this line !
  for i in 1:2
	   x = 4
  end
  return x
end
