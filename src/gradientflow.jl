function flowode!(dz,z,p,t)
    x = z[1]
    y = z[2]

    dz[1] = -(4*x^3+4*x*y-44*x+2*x+2*y^2-14)
    dz[2] = -(4*y*x+4*y^3-28*y+2*x^2+y-22)
end

using OrdinaryDiffEq

Tf = 200.0;
tspan = (0.0, Tf)
x0 = [-1,-1] # zeros(2) # [0.5,0.1]
pars = [-7,-11]
alg = Trapezoid() # Tsit5() # KenCarp5()
tsave = 0.1;
prob = ODEProblem(flowode!,x0, tspan)
sol = solve(prob, alg, p=pars, save_everystep=false)

data = Array(sol)'

using Plots
plot(data[:,1], data[:,2])

data[end,:]

using LinearAlgebra

xs = -5:0.1:5
ys = -5:0.1:5

cloud = zeros(length(xs),length(ys))

p_min = [3.0 2.0 1;
        -2.805118 3.131312 2;
        -3.779310 -3.283186 3;
        3.584428 -1.848126 4];


for (yi,y) in enumerate(ys), (xi,x) in enumerate(xs)
    x0 = [x,y] 
    prob = ODEProblem(flowode!,x0, tspan)
    sol = solve(prob, alg, p=pars, save_everystep=false)

    data = sol[end]

    for i = 1 : 4
        e = data-p_min[i,1:2]    
        if norm(e) < 0.1
            cloud[xi,yi] = p_min[i,3]
        end
    end
end

heatmap(cloud)

