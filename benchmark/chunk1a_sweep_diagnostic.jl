# Error/refinement sweeps against the analytic oracle (diagnostic, not a test).
using Printf, ForwardDiff, LinearAlgebra
using SciMLBase: ODEProblem, solve
using ADTypes: AutoFiniteDiff
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: QNDF
include(joinpath(@__DIR__, "..", "test", "chunk1a_helpers.jl"))
Jan = ForwardDiff.jacobian(analytic_g, THETA_NOMINAL)
colscale = [maximum(abs.(Jan[:, j])) for j in 1:5]
names = ["a", "c", "b", "y10", "y20"]
function report(label, J)
    D = abs.(J .- Jan)
    scaled = D ./ colscale'
    idx = argmax(scaled)
    @printf("  %-34s max|dJ|=%.3e  max col-scaled=%.3e  worst (row %d, col %s): J=%.6e ref=%.6e\n",
        label, maximum(D), scaled[idx], idx[1], names[idx[2]], J[idx], Jan[idx])
    println("      col-scaled per column: ", join([@sprintf("%s=%.2e", names[j], maximum(scaled[:, j])) for j in 1:5], " "))
end
println("analytic Jacobian (rows: y1,y2 @ early; y1,y2 @ interp; y1,y2 @ late):"); display(Jan); println()
println("column scales: ", colscale)
for (nm, alg) in (("Rodas5P", Rodas5P()), ("Rodas5P(FiniteDiff)", Rodas5P(autodiff = AutoFiniteDiff())), ("QNDF", QNDF()))
    println("== $nm")
    for (at, rt) in ((1e-6, 1e-5), (1e-8, 1e-7), (1e-10, 1e-9), (1e-12, 1e-11), (1e-13, 1e-12))
        g = th -> numeric_g(th, alg; abstol = at, reltol = rt)
        println(" tol abstol=$at reltol=$rt")
        pr = maximum(abs.(g(THETA_NOMINAL) .- analytic_g(THETA_NOMINAL)))
        @printf("  primal max abs err vs analytic = %.3e\n", pr)
        Jfd = ForwardDiff.jacobian(g, THETA_NOMINAL)
        report("ForwardDiff(solve) vs analytic", Jfd)
        for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
            report("central FD h=$h vs analytic", fd_jacobian(g, THETA_NOMINAL, h))
        end
    end
end
