using Random
using Statistics
using GLMakie
# function fractalsurface(n, H, seed)



function fractalsurface(n, H, seed)
    rand_gen = Xoshiro(seed)

    N = 2^n

    L = zeros(N+1,N+1)

    R = randn(rand_gen, 3,3)

    initial = [1,2^(n-1),2^n+1]

    for i in 1:3
        for j in 1:3
            L[initial[i], initial[j]]=R[i,j]
        end
    end

    idx = findall(!iszero, L)
    r = Int(sqrt(length(idx)))


    for l=1:n-1
        mid = zeros(Int, r-1)
        for i=1:r-1
            mid[i]=Int((idx[i][1]+idx[i+1][1])÷2);
            for j=1:r
                L[mid[i],idx[j][1]] = mean([L[idx[i][1],idx[j][1]], L[idx[i+1][1], idx[j][1]]]) + randn(rand_gen)*2^(-H*(l+1));
                L[idx[j][1], mid[i]] = mean([L[idx[j][1], idx[i][1]], L[idx[j][1], idx[i+1][1]]]) + randn(rand_gen)*2^(-H*(l+1));
            end
        end
        for i=1:r-1
            for j=1:r-1
                L[mid[i],mid[j]] = mean([L[idx[i][1],idx[j][1]], L[idx[i][1], idx[j+1][1]], L[idx[i+1][1], idx[j][1]], L[idx[i+1][1], idx[j+1][1]]])+randn(rand_gen)*2^(-H*(l+1));
            end
        end

        idx = findall(!iszero, L)

        r = Int(sqrt(length(idx)));

    end

    R = zeros(r,r)

    for i=1:r
        for j=1:r
            R[i,j] = L[idx[i][1],idx[j][1]];
        end
    end

    return R
end


function fractalasperity(n, H)

    seed = rand(Int16)

    surf = fractalsurface(n, H, seed)
    fractal = surf .< -1
    return fractal, seed
end


function fractalasperity(n, H, seed)

    surf = fractalsurface(n, H, seed)
    fractal = surf .< -1
    return fractal
end

# seed = 2
n = 9
H = 0.45

fig = Figure()

ax = Axis(fig[1, 1], aspect=DataAspect())

sg = SliderGrid(
    fig[2, 1],
    (label = "n", range = 2:1:12, startvalue = 5),
    (label = "H", range = 0:0.01:1, format = "{:.2f}", startvalue = 0.45),
    (label = "Seed", range = -100:1:100, startvalue = 0),
    width = 350,
    tellheight = false)

sliderobservables = [s.value for s in sg.sliders]


fractal = lift(sliderobservables...) do slvalues...
    fractalasperity(slvalues...)


end

heatmap!(ax, fractal)


display(fig)
