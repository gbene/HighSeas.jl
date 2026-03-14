using CUDA
using LinearAlgebra
using BenchmarkTools
using FastBroadcast

function foo(temp, x, mask, total_nodes)
    @.. thread=true temp = x*mask
    return sum(temp)/total_nodes
end

function foo(x, mask, total_nodes)
    return dot(x, mask) / total_nodes
end

function foo_section(res, temp, x, mask)
    @.. thread=true temp = x*mask
    res[5,:] =  sum(temp, dims=1)
    return res
end

function foo_section(res, x, mask)
    res[5,:] .= vec(mapreduce(*, +, x, mask; dims=1))
    return res
end

N = 1024

temp = zeros(N, N)
x = rand(N, N)
mask = @. Int8(x >= 0.3)
section_mask = zeros(Int8, N, N)
section_mask[5, :] .= 1

temp_cu = CUDA.CuMatrix{Float64}(temp)
x_cu = CUDA.CuMatrix{Float64}(x)
mask_cu = CUDA.CuMatrix{Int8}(mask)
section_mask_cu = CUDA.CuMatrix{Int8}(section_mask)


total_nodes = sum(mask)

res = zeros(N,N)
res_cu = CUDA.CuMatrix{Float64}(res)

# @btime foo($temp, $x, $mask, $total_nodes)
# @btime CUDA.@sync foo($temp_cu, $x_cu, $mask_cu, $total_nodes)
# @btime foo($x, $mask, $total_nodes)
# @btime CUDA.@sync foo($x_cu, $mask_cu, $total_nodes)

@btime foo_section($res, $temp, $x, $section_mask)
@btime CUDA.@sync foo_section($res_cu, $temp_cu, $x_cu, $section_mask_cu)
@btime foo_section($res, $x, $section_mask)
@btime CUDA.@sync foo_section($res_cu, $x_cu, $section_mask_cu)
