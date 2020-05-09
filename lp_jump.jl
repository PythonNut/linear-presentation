using JuMP, Gurobi, Combinatorics

function build_paths(MOD, n, bound, m)
    @variable(MOD, 0 <= paths[1:n, 1:m] <= bound, Int)
    @variable(MOD, inds[1:n, 1:m], Bin)
    @variable(MOD, ne_inds[1:n, 1:m], Bin)

    for i in 1:n
        for j in 1:m
            @constraint(MOD, bound*inds[i, j] >= paths[i, j])
            @constraint(MOD, paths[i, j] >= inds[i, j])
        end

        for j in 1:m-1
            @constraint(MOD, inds[i, j] >= inds[i, j+1])

            hot = -2*bound*(2 - inds[i, j] - inds[i, j+1])
            @constraint(MOD, paths[i, j] - paths[i, j+1] + 2*bound*ne_inds[i, j] >= 1 + hot)
            @constraint(MOD, paths[i, j+1] - paths[i, j] + 2*bound*(1-ne_inds[i, j]) >= 1 + hot)
        end
    end

    @variable(MOD, parity[1:n], Bin)
    return paths, parity, inds
end

function require_planarity(MOD, paths, parity, inds, n, bound, m)
    @variable(MOD, pln_inds[1:n, 1:m, 1:n, 1:m, 1:2], Bin)
    @variable(MOD, 0 <= maxs[1:n, 1:m] <= bound, Int)
    @variable(MOD, max_inds[1:n, 1:m], Bin)

    for i in 1:n
        for j in 1:m-1
            @constraint(MOD, maxs[i, j] >= paths[i, j])
            @constraint(MOD, maxs[i, j] >= paths[i, j+1])

            @constraint(MOD, paths[i, j] + 2*bound*max_inds[i, j] >= maxs[i, j])
            @constraint(MOD, paths[i, j+1] + 2*bound*(1-max_inds[i, j]) >= maxs[i, j])
        end
    end

    for (i, k) in combinations(1:n, 2)
        for j in 1:m-1, l in 1:m-1
            b = maxs[i, j]
            a = paths[i, j] + paths[i, j+1] - b
            d = maxs[k, l]
            c = paths[k, l] + paths[k, l+1] - d
            M = 2*bound + 1

            flipi = if j%2 == 1 parity[i] else (1 - parity[i]) end
            flipk = if l%2 == 1 parity[k] else (1 - parity[k]) end

            hot = M*(4 - inds[i, j] - inds[i, j+1] - inds[k, l] - inds[k, l+1])
            both_bot = M * (flipi + flipk)
            both_top = M * (2 - flipi - flipk)

            inside = pln_inds[i,j,k,l,1]
            branch = pln_inds[i,j,k,l,2]

            @constraint(MOD, b <= c-1 + M*(1-inside) + M*branch + hot + both_top)
            @constraint(MOD, b <= c-1 + M*(1-inside) + M*branch + hot + both_bot)

            # c < d < a < b
            @constraint(MOD, d <= a-1 + M*(1-inside) + M*(1-branch) + hot + both_top)
            @constraint(MOD, d <= a-1 + M*(1-inside) + M*(1-branch) + hot + both_bot)

            # a < c < d < b
            @constraint(MOD, a <= c-1 + M*inside + M*branch + hot + both_top)
            @constraint(MOD, d <= b-1 + M*inside + M*branch + hot + both_top)
            @constraint(MOD, a <= c-1 + M*inside + M*branch + hot + both_bot)
            @constraint(MOD, d <= b-1 + M*inside + M*branch + hot + both_bot)

            # c < a < b < d
            @constraint(MOD, c <= a-1 + M*inside + M*(1-branch) + hot + both_top)
            @constraint(MOD, b <= d-1 + M*inside + M*(1-branch) + hot + both_top)
            @constraint(MOD, c <= a-1 + M*inside + M*(1-branch) + hot + both_bot)
            @constraint(MOD, b <= d-1 + M*inside + M*(1-branch) + hot + both_bot)
        end
    end
    return pln_inds
end

function require_connections(MOD, semiarcs, paths, parity, inds, step_size, shift, m, bound)
    for (i, (a, da, b, db)) in enumerate(semiarcs)
        @constraint(MOD, paths[i, m] == 0)

        # Taking care of the "a" conditions is easy
        a_pos = step_size * a + shift
        if da == 0
            @constraint(MOD, parity[i] == 1)
        elseif da == 2
            @constraint(MOD, parity[i] == 0)
        elseif da == 1
            a_pos += 1
            if db == 3
                @constraint(MOD, parity[i] == 1)
            end
        end

        @constraint(MOD, paths[i, 1] == a_pos)

        # Now, for the "b" conditions
        b_pos = step_size * b + shift
        if db == 3
            b_pos -= 1
        end

        for j in 1:m-1
            this_slice = inds[i, j]
            next_slice = inds[i, j+1]
            @constraint(MOD, 2*bound * ((1 - this_slice) + next_slice) >= paths[i, j] - b_pos)
            @constraint(MOD, 2*bound * ((1 - this_slice) + next_slice) >= b_pos - paths[i, j])

            flip = if j%2 == 1 parity[i] else (1 - parity[i]) end
            if db == 0
                @constraint(MOD, 1-flip >= this_slice - next_slice)
            elseif db == 2
                @constraint(MOD, flip >= this_slice - next_slice)
            end
        end
    end
end


semiarcs = [
    (1, 1, 2, 3), (2, 1, 3, 3), (3, 1, 4, 3), (4, 1, 5, 3),
    (5, 1, 6, 3), (6, 1, 7, 3), (7, 1, 8, 3), (8, 1, 9, 3),
    (9, 1, 1, 0), (1, 2, 2, 2), (2, 0, 10, 3), (10, 1, 11, 3),
    (11, 1, 12, 3), (12, 1, 6, 0), (6, 2, 5, 2), (5, 0, 12, 2),
    (12, 0, 13, 3), (13, 1, 10, 0), (10, 2, 3, 0), (3, 2, 14, 3),
    (14, 1, 15, 3), (15, 1, 16, 3), (16, 1, 9, 2), (9, 0, 8, 0),
    (8, 2, 7, 2), (7, 0, 13, 2), (13, 0, 11, 0), (11, 2, 4, 0),
    (4, 2, 16, 0), (16, 2, 15, 2), (15, 0, 14, 0), (14, 2, 1, 3)
]

n = 16
m = 4

MOD = Model(Gurobi.Optimizer)
set_optimizer_attribute(MOD, "Seed", 3)
set_optimizer_attribute(MOD, "MIPFocus", 3)
set_optimizer_attribute(MOD, "Presolve", 2)

paths, parity, inds = build_paths(MOD, length(semiarcs), n^2, m)
pln_inds = require_planarity(MOD, paths, parity, inds, length(semiarcs), n^2 + 2, m)
require_connections(MOD, semiarcs, paths, parity, inds, n, -n+2, m, n^2 + 2)

@objective(MOD, MOI.MIN_SENSE, sum(inds))
