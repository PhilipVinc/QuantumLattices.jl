using Documenter, NeuralQuantum

makedocs(
    modules   = [NeuralQuantum],
    format    = :html,
    sitename  = "Quantumoptics.jl",
    authors   = "Filippo Vicentini",
    pages     = [
            "Home"          => "index.md"
    ]
)

deploydocs(
    repo   = "github.com/PhilipVinc/QuantumOptics.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
