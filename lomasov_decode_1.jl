using Random
#умножение по модулю 2
function mul_mod2(A, B)

    if ndims(A) == 2 && ndims(B) == 1
        m, n = size(A)
        res = zeros(Int, m)
        for i in 1:m
            s = 0
            for j in 1:n
                s += A[i, j] * B[j]
            end
            res[i] = s % 2
        end
        return res
    end
  
    if ndims(A) == 1 && ndims(B) == 2
        n, p = size(B)
        res = zeros(Int, p)
        for j in 1:p
            s = 0
            for i in 1:n
                s += A[i] * B[i, j]
            end
            res[j] = s % 2
        end
        return res
    end

    error("error")
end
function inv_mod2(A)
    n = size(A, 1)

    I = zeros(Int, n, n)#ед. матрица
    for i in 1:n
        I[i, i] = 1
    end

    A_ext = hcat(copy(A), I)

    for i in 1:n
        if A_ext[i, i] == 0
            pivot = findfirst(x -> x == 1, A_ext[i+1:end, i])
            if pivot === nothing
                return Matrix{Int}(undef, 0, 0), false
            end
            pivot += i
            tmp = copy(A_ext[i, :])
            A_ext[i, :] = A_ext[pivot, :]
            A_ext[pivot, :] = tmp
        end

        for j in 1:n
            if j != i && A_ext[j, i] == 1
                for k in 1:2n
                    A_ext[j, k] = (A_ext[j, k] + A_ext[i, k]) % 2
                end
            end
        end
    end

    return A_ext[:, n+1:end], true
end

#decodeISD
function decodeISD(G, y, t; niter = -1)
    k, n = size(G)
    iter = 0

    while niter < 0 || iter < niter
        iter += 1

        I = sort(rand(1:n, k))

        GI = G[:, I]
        GI_inv, ok = inv_mod2(GI)
        if !ok
            continue
        end

        m_candidate = mul_mod2(transpose(GI_inv), y[I])
        y_hat = mul_mod2(m_candidate, G)

        e = [(y[i] + y_hat[i]) % 2 for i in 1:n]

        if sum(e) == t
            println("решение найдено,на итерации $iter")
            return m_candidate, e
        end
    end

    println("решение не найдено, $iter итераций")
    return nothing
end

println("щапуск теста\n")

G = [
    1 0 0 1 1 0;
    0 1 0 0 1 1;
    0 0 1 1 0 1
]

y = [1, 1, 1, 0, 1, 0]

t = 2



#запуск
result = decodeISD(G, y, t; niter=500)
println("\nРезультат:")
println(result)
