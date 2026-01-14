using Random

function get_combos(n, p)
    combos = Vector{Vector{Int}}()
    function backtrack(idx, w, curr)
        if w == 0
            push!(combos, copy(curr))
            return
        end
        for i in idx:n
            curr[i] = 1
            backtrack(i + 1, w - 1, curr)
            curr[i] = 0
        end
    end
    backtrack(1, p, zeros(Int, n))
    return combos
end

function decode_stern(H_in, t_target, p_param, l_param, niter=-1)
    r, n = size(H_in)
    k = n - r
    k1 = k ÷ 2
    k2 = k - k1
    iter_count = 0

    while niter < 0 || iter_count < niter
        iter_count += 1
        perm = randperm(n)
        
        M = zeros(Int, r, n)
        for i in 1:r
            for j in 1:n
                M[i, j] = H_in[i, perm[j]]
            end
        end

        pivot_cols = Int[]
        row = 1
        for j in 1:n
            row > r && break
            sel = -1
            for i in row:r
                if M[i, j] == 1
                    sel = i
                    break
                end
            end
            if sel != -1
                M[row, :], M[sel, :] = M[sel, :], M[row, :]
                for i in 1:r
                    if i != row && M[i, j] == 1
                        for col in 1:n
                            M[i, col] = M[i, col] != M[row, col] ? 1 : 0
                        end
                    end
                end
                push!(pivot_cols, j)
                row += 1
            end
        end

        length(pivot_cols) < r && continue

        other_cols = Int[]
        for j in 1:n
            if !(j in pivot_cols) push!(other_cols, j) end
        end

        P = zeros(Int, r, k)
        for i in 1:r
            for j in 1:k
                P[i, j] = M[i, other_cols[j]]
            end
        end

        table = Dict{Vector{Int}, Vector{Vector{Int}}}()
        for e1 in get_combos(k1, p_param)
            val = zeros(Int, r)
            for j in 1:k1
                if e1[j] == 1
                    for i in 1:r
                        val[i] = val[i] != P[i, j] ? 1 : 0
                    end
                end
            end
            key = val[1:l_param]
            if !haskey(table, key) table[key] = [] end
            push!(table[key], e1)
        end

        for e2 in get_combos(k2, p_param)
            val = zeros(Int, r)
            for j in 1:k2
                if e2[j] == 1
                    for i in 1:r
                        val[i] = val[i] != P[i, k1+j] ? 1 : 0
                    end
                end
            end
            
            key = val[1:l_param]
            if haskey(table, key)
                for e1 in table[key]
                    e_syst = zeros(Int, r)
                    for i in 1:r
                        v1_i = 0
                        for j in 1:k1
                            if e1[j] == 1 v1_i = v1_i != P[i, j] ? 1 : 0 end
                        end
                        e_syst[i] = v1_i != val[i] ? 1 : 0
                    end

                    if sum(e_syst) + sum(e1) + sum(e2) == t_target
                        res_perm = vcat(e_syst, e1, e2)
                        e_final = zeros(Int, n)
                        final_order = vcat(pivot_cols, other_cols)
                        for i in 1:n
                            e_final[perm[final_order[i]]] = res_perm[i]
                        end
                        return e_final, iter_count
                    end
                end
            end
        end
    end
    return nothing, iter_count
end

# --- ВХОДНЫЕ ДАННЫЕ ---

H = [
    1 1 1 1 1 1 1 1;
    0 0 0 0 1 1 1 1;
    0 0 1 1 0 0 1 1;
    0 1 0 1 0 1 0 1
]

t=4
p=1
l=1
niter=500
ans, iterations = decode_stern(H, t, p, l,niter)
println("Количество итераций: ", iterations)
println("Вектор: ", ans)