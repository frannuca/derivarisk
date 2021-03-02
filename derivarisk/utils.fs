namespace DerivaRisk

module Utils=
    open System
    open FSharp.Stats


    let randu N=
        let u = fun _ -> Distributions.Continuous.Uniform.Sample 0.0 1.0
        Array.init N u

    let randn N=
        let u = fun _ -> Distributions.Continuous.Normal.Sample 0.0 1.0
        Array.init N u