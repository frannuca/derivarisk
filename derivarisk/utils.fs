namespace DerivaRisk

module Utils=
    open System
    open MathNet.Numerics.LinearAlgebra
    open FSharp.Stats
    open Deedle
    open Deedle.Math


    let randu N=
        let u = fun _ -> Distributions.Continuous.Uniform.Sample 0.0 1.0
        Array.init N u

    let randn N=
        let u = fun _ -> Distributions.Continuous.Normal.Sample 0.0 1.0
        Array.init N u

    let frame2Matrix  rows columns frame=
        frame |> Frame.sliceCols columns |> Frame.sliceRows rows |> Deedle.Math.Matrix.ofFrame

    let matrix2Frame rows columns (m:MathNet.Numerics.LinearAlgebra.Matrix<float>) =
       Frame.ofMatrix rows columns m

    let writeCsv(mcpaths:float[][][],filepath)=        
        use stream = new System.IO.StreamWriter(filepath,false)
        mcpaths |> Array.iteri(fun nsim assets -> assets|>Array.iteri(fun nasset asset -> stream.WriteLine(nsim.ToString()+","+nasset.ToString()+","+String.Join(",",asset))))