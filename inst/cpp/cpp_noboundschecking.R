# gsub(".at\\((.*?)\\)", "[\\1]", "abab.at(kk) = f(x)")
# gsub(".at\\((.*?)\\)", "[\\1]", "if(UU[i]>0)")

cpp1 <- readLines("./src/gfilmm.cpp")
cpp2 <- gsub(".at\\((.*?)\\)", "[\\1]", cpp1)
cpp2 <- gsub(".at\\((.*?)\\)", "[\\1]", cpp2)
# not perfect yet

writeLines(cpp2, "./inst/cpp/gfilmm_NBC.cpp")
