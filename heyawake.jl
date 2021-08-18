using JuMP, CPLEX, SparseArrays

function readFile(sfile)
   sf = open(sfile,"r")

   s = read(sf,String)

   close(sf)

   s = parse.(Int,split(s))

   n = s[1]  # dimensao do tabuleiro  
   r = s[2]  # nro de regioes (salas) 

   k = 3;
   rooms = zeros(Int,r,5)
   for i = 1 : r
         rooms[i,1:5] =  s[k:k+4] 
         k += 5
   end

   return n, rooms

end

# =========================

function buildModel(n,rooms)

   model = Model(with_optimizer(CPLEX.Optimizer, CPX_PARAM_TILIM=1800))  #CPX_PARAM_SCRIND=false, 

   @variable(model, x[i = 1:n, j = 1:n], Bin)  

   @objective(model,Min,sum(x[i,j] for i = 1:n for j = 1 : n))

# =====================================  Restricao Vizinhos Ortogonais ===========================================
   for i = 1 : n
      for j = 1 : n
         if( j + 1 <= n ) 
            @constraint(model,x[i,j] + x[i,j+1] <= 1)
         end
         if( i + 1 <= n ) 
            @constraint(model,x[i,j] + x[i+1,j] <= 1)
         end
      end
   end

# =====================================  Restricao Nro Obrigatorio de Celulas Sombreadas por Regioes ===============
   r = length(rooms[:,1])
   for k = 1 : r
      if rooms[k,5] > -1
          @constraint(model, sum(  x[i,j] for i = rooms[k,1] : rooms[k,3] for j = rooms[k,2] : rooms[k,4] )  == rooms[k,5] )
      end                    
   end
# =====================================  Restricao Linha Reta Cruzando Tres Regioes ===========================================
   tableau = zeros(Int,n,n)
   for k = 1 : r
      for i = rooms[k,1] : rooms[k,3]
         for j = rooms[k,2] : rooms[k,4]
            tableau[i,j] = k 
         end
      end
   end

   for i = 1 : n             #  considera linhas
      cr = tableau[i,1]
      inic  = 0 
      k = 1 
      j = 2
      while  j <= n 
          if cr != tableau[i,j]
             if k == 1
                inic =  j-1
             end
             cr = tableau[i,j]
             k += 1
          end
          if k == 3
             @constraint(model, sum(x[i,jl]  for jl = inic : j)  >= 1 )
             inic += 1
             j = inic +1
             cr = tableau[i,inic]
             k = 1 
          else  
             j = j + 1 
          end
      end
   end

   for j = 1 : n      # considera colunas
      cr = tableau[1,j]
      inic  = 0 
      k = 1 
      i = 2
      while  i <= n 
          if cr != tableau[i,j]
             if k == 1
                inic =  i-1
             end
             cr = tableau[i,j]
             k += 1
          end
          if k == 3
             @constraint(model, sum(x[il,j]  for il = inic : i)  >= 1 )
             inic += 1
             i = inic +1
             cr = tableau[inic,j]
             k = 1 
          else  
             i = i + 1 
          end
      end
   end


# =====================================  Restricao continuidade area com celulas nao sombreadas =======================================
   @variable(model, yR[i = 1:n,   j = 1:n-1, i, j+1] >=0)  
   @variable(model, yD[i = 1:n-1, j = 1:n  , i+1, j] >=0)  
   @variable(model, yL[i = 1:n,   j = 2:n  , i  , j-1] >=0)  
   @variable(model, yU[i = 2:n,   j = 1:n  , i-1, j] >=0 )   
   @variable(model, f[i = 1:2], Bin)  
   
   @constraint(model,f[1] + f[2] <= 2 - x[1,1])
   @constraint(model,yD[1,1,2,1] <= f[1] * n^2 )
   @constraint(model,yR[1,1,1,2] <= f[2] * n^2 )

   @constraint(model,yD[1,1,2,1] + yR[1,1,1,2] == n^2 - 1 - sum(x[i,j] for i = 1 : n for j = 1 : n) + x[1,1])
   @constraint(model,yR[1,n-1,1,n] - yL[1,n,1,n-1] + yU[2,n,1,n] - yD[1,n,2,n] == 1 - x[1,n])
   @constraint(model,yL[n,2,n,1] - yR[n,1,n,2] + yD[n-1,1,n,1] - yU[n,1,n-1,1] == 1 - x[n,1])
   @constraint(model,yR[n,n-1,n,n] - yL[n,n,n,n-1] + yD[n-1,n,n,n] - yU[n,n,n-1,n] == 1 - x[n,n])

   for j = 2 : n-1
      @constraint(model,yR[1,j-1,1,j] + yU[2,j,1,j] + yL[1,j+1,1,j] - yD[1,j,2,j] - yL[1,j,1,j-1] - yR[1,j,1,j+1] == 1 - x[1,j])
      @constraint(model,yR[n,j-1,n,j] + yD[n-1,j,n,j] + yL[n,j+1,n,j] - yU[n,j,n-1,j] - yL[n,j,n,j-1] - yR[n,j,n,j+1] == 1 - x[n,j])
   end
   for i = 2 : n-1
      @constraint(model,yL[i,2,i,1] + yD[i-1,1,i,1] + yU[i+1,1,i,1] - yR[i,1,i,2] - yU[i,1,i-1,1] - yD[i,1,i+1,1] == 1 - x[i,1])
      @constraint(model,yR[i,n-1,i,n] + yD[i-1,n,i,n] + yU[i+1,n,i,n] - yL[i,n,i,n-1] - yU[i,n,i-1,n] - yD[i,n,i+1,n] == 1 - x[i,n])
   end

   for i = 2 : n-1 
      for j = 2 : n-1
         @constraint(model, yD[i-1,j,i,j] + yU[i+1,j,i,j] + yL[i,j+1,i,j] + yR[i,j-1,i,j]   
                          - yU[i,j,i-1,j] - yD[i,j,i+1,j] - yR[i,j,i,j+1] - yL[i,j,i,j-1]  == 1 - x[i,j])
      end
   end

   @constraint(model,yR[1,n-1,1,n] + yU[2,n,1,n]  <= n^2 * (1 - x[1,n]))
   @constraint(model,yL[n,2,n,1] + yD[n-1,1,n,1]  <= n^2 * (1 - x[n,1]))
   @constraint(model,yR[n,n-1,n,n] + yD[n-1,n,n,n] <= n^2 * (1 - x[n,n]))

   for j = 2 : n-1
      @constraint(model,yR[1,j-1,1,j] + yU[2,j,1,j] + yL[1,j+1,1,j]  <= n^2 * (1 - x[1,j]))
      @constraint(model,yR[n,j-1,n,j] + yD[n-1,j,n,j] + yL[n,j+1,n,j] <= n^2 * (1 - x[n,j]))
   end
   for i = 2 : n-1
      @constraint(model,yL[i,2,i,1] + yD[i-1,1,i,1] + yU[i+1,1,i,1]  <= n^2 * (1 - x[i,1]))
      @constraint(model,yR[i,n-1,i,n] + yD[i-1,n,i,n] + yU[i+1,n,i,n] <= n^2 * (1 - x[i,n]))
   end

   for i = 2 : n-1 
      for j = 2 : n-1
         @constraint(model, yD[i-1,j,i,j] + yU[i+1,j,i,j] + yL[i,j+1,i,j] + yR[i,j-1,i,j] <= n^2 * (1 - x[i,j]))
      end
   end



   JuMP.optimize!(model)


   tableau = zeros(n,n);
   for i = 1 : n
      for j =  1 : n
         if JuMP.value(x[i,j]) > 0.5 
            println("x[$i,$j] = ",JuMP.value(x[i,j]))
            tableau[i,j] = 1
         end
      end
   end
   
   conta  = checkSolu(tableau)
   if conta == 0
      println("Solucao ok")
   else
      println("$conta celulas ilhadas")
   end

   

#   println("y[1,1,1,2] = ",JuMP.value(yR[1,1,1,2]))

   open("model.lp", "w") do f
      print(f, model)
   end

end

# =========================
function checkSolu(tableau)
   n = length(tableau[:,1])
   i = 1
   j = 1

   while tableau[i,j] > 0 
      j += 1
      if j > n 
         i += 1
         j = 1
      end
   end   

   verificaViz(tableau,n,i,j)

   ok = 0  
   for i = 1 : n 
      for j = 1 : n
         if tableau[i,j] == 0
            ok += 1
         end
      end
   end

   return ok
end

# =========================
function verificaViz(tableau,n,i,j)
   if tableau[i,j] == 0
      tableau[i,j] = 1
      if  i - 1 > 0 
         verificaViz(tableau,n,i-1,j)
      end
      if  i + 1 <= n 
         verificaViz(tableau,n,i+1,j)
      end
      if  j - 1 > 0 
         verificaViz(tableau,n,i,j-1)
      end
      if  j + 1 <= n 
         verificaViz(tableau,n,i,j+1)
      end
   end
end

# =========================

n, rooms = readFile(ARGS[1])

buildModel(n,rooms)


