namespace DerivaRisk
open Deedle
module MCSimCube=
    open MathNet.Numerics.LinearAlgebra
    open MathNet.Numerics.Random
   

    (*
        <summary>
            Generation of Monte Carlo paths.
        </summary>
        <param name="rho"> correlation matrix amongh assets</param>
        <param name="npaths">number of simulation paths</param>
        <param name="timeSteps">number of time steps included in each simulation path</param>
        <param name="number_assets">number of correlated assets included in each time step</param>
        <returns>3D array with the following axis definitions:
                    -[0]-> simulation path index
                    -[1]-> asset index
                    -[2]-> time step index
        </returns>
    *)
    let generate_cube(rho:Matrix<float>)(seed:int)(npaths:int)(timeSteps:int)(number_assets:int):float[][][]=
        let rndgen = new MersenneTwister(seed,true)
        let norm = new MathNet.Numerics.Distributions.Normal(rndgen)
        let C = rho.Cholesky().Factor
        //printfn "Cholesky Factor %A" C
        let path_generator _ =
            let mcpath = [|0 .. number_assets-1|] |> Array.map(fun _ -> let aux = Array.zeroCreate(timeSteps)
                                                                        norm.Samples(aux) |> ignore
                                                                        aux                                                         
                                                )   |> Matrix<float>.Build.DenseOfRowArrays

            (C*mcpath).ToRowArrays()

        let simcube = [|0 .. npaths-1|] |> Array.map(path_generator)
        simcube
        

    let apply_payoff(cube:float[][][])(payoff:float[][]->float):float[]=
        cube
        |>Array.map(payoff)
       