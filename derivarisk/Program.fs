// Learn more about F# at http://fsharp.org

open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
[<EntryPoint>]
let main argv =

    let nsim=10000
    let ntime=100
    let T=0.25
    let K=90.0
    let hestonpar = {   HestonParams.dividends=0.02;
                        rate=0.03;
                        kappa=6.2;
                        theta=0.06;
                        sigma=0.5;
                        rho= -0.7;
                        lambda=0.0;
                        S0=100.0;
                        V0=0.03;
                        dt=T/float ntime;
                        gamma_1=0.5;
                        gamma_2=0.5;
                        threshold_phi=1.5
                     }

    let heston_call = Heston.heston_analytical(hestonpar.S0,
                                               K,
                                               hestonpar.kappa,
                                               0.0,
                                               T,
                                               hestonpar.rate,
                                               hestonpar.dividends,
                                               hestonpar.rho,
                                               hestonpar.sigma,
                                               hestonpar.theta,
                                               hestonpar.V0,
                                               true)

    let hestonparams=Map([Heston.AssetName("A"),hestonpar;
                          Heston.AssetName("B"),{hestonpar with rho=0.7};
                          Heston.AssetName("C"),{hestonpar with sigma=0.15}])


    let rs= 0.5
                          
    let assets = [Heston.AssetName "A";Heston.AssetName "B";Heston.AssetName "C"]
                          
    let rho = Matrix<float>.Build.DenseOfRowArrays(
                                        [|
                                            [|1.0; rs;rs|];
                                            [|rs; 1.0;rs|];
                                            [|rs; rs;1.0|]
                                        |])|> Utils.matrix2Frame assets assets

    let rv = rs
    let rho_vol = Matrix<float>.Build.DenseOfRowArrays(
                        [|
                            [|1.0; rv;rv|];
                            [|rv; 1.0;rv|];
                            [|rv; rv;1.0|]
                        |])|> Utils.matrix2Frame assets assets
        
    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let mcpaths, cube = duration (fun () -> Heston.computeMCPaths(rho,rho_vol,nsim,ntime,hestonparams))
    Utils.writeCsv(mcpaths,"/Users/fran/data/mcpath_multiple_asset_different.csv")

    let spath = mcpaths.[0]
    let ld = spath
            |> Array.map(fun path ->[|1 .. path.Length-1|]
                                    |> Array.map(fun n -> path.[n]/path.[n-1]-1.0))
                            

    let rhoS = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(ld)
    printfn "Correlation rhoS \n %A" rhoS

    let rhosCube = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(cube.[0])
    printfn "Correlation rhoCube \n %A" rhosCube

    printfn "Difference %f" ((rhoS-rhosCube).ToColumnMajorArray() |> Array.maxBy(Math.Abs))
    0 // return an integer exit code
