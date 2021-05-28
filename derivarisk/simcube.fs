namespace DerivaRisk.model
open InstrumentData
open DerivaRisk.data
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
module MonteCarloCube=
    let generate_cube(config:SimulationCube)(seed:int):string[]*float[][][]=
        let assets,rho = config.assets_correlation |> Extractors.messageMatrix2Matrix        
        let number_assets = int(config.number_of_assets)
        let timeSteps = int(config.number_of_time_steps)
        let npaths = int(config.number_of_simulations)

        let rndgen = new MersenneTwister(int(seed),true)
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
        assets,simcube