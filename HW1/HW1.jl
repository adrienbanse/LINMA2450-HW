using JuMP, Gurobi, XLSX

xf = XLSX.readxlsx("toy-stores.xlsx")

data = xf["Sheet1"]

model = Model(with_optimizer(Gurobi.Optimizer))
M = 100
N = 24

set_optimizer_attribute(model, "Presolve", 0)
set_optimizer_attribute(model, "Heuristics", 0)
set_optimizer_attribute(model, "Cuts", 0)

# Variables
@variable(model, 0 <= x[1:M,1:N]) # + Constraint (4)
@variable(model, 0 <= y[1:N],Bin) # + Constraint (5)

# Compute costs C_ij and fixed costs f_j
C = zeros(Float64, (M,N))
f = zeros(Float64,N)
for j in 1:N
    lat2 = deg2rad(data[j+2,11])
    long2 = deg2rad(data[j+2,10])
    f[j]=data[j+2,12]
    for i in 1:M
        lat = deg2rad(data[i+2,4]) 
        long = deg2rad(data[i+2,3])
        if lat==lat2 && long ==long2 
            C[i,j]=0
        else
            C[i,j]=2*asin(sqrt( sin((lat2-lat)/2)^2 + cos(lat)*cos(lat2)*sin((long2-long)/2)^2 ))*3958.8 *data[i+2,5]
            #C[i,j]= acos(sin(lat)*sin(lat2)+cos(lat)*cos(lat2)*cos(long-long2))
        end
    end
end

# Objective function
@objective(model, Min, sum(C .* x) + sum(f .* y))

# Constraint (2)
@constraint(model, sum(x, dims=2) .== ones(Float64, M))

# Constraint (3)
@constraint(model, [i = 1:M], x[i,:] .<= y)

# Constraint (6)
#@constraint(model, sum(x, dims=1) .<= M * transpose(y))

println(optimize!(model))

#findall(x->x!=-0.0, JuMP.value.(y))
#findall(y->(y!=1.0 && y!=0.0), JuMP.value.(x))
#xsum = sum(JuMP.value.(x),dims=1)
println(xsum[findall(x->x!=0.0,xsum)])
#-> A savoir : Buffalo, Houston, Minneapolis, Oakland, Orlando, Paradise CDP, Winston-Salem
