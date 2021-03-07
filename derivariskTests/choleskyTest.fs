module HestonTests
open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open NUnit.Framework
open System.IO

[<SetUp>]
let Setup () =
    ()

[<Test>]
let TestCholesky () =
    let rho = Matrix<float>.Build.DenseOfRowArrays([|
        [|1.0; 0.75|];
        [|0.75; 1.0|]
       |])

    let nsim=2
    let ntime=100
    let nassets=2
    let seed=42
    let cube = MCSimCube.generate_cube  rho seed nsim ntime nassets
    

    let crho = ([for path in cube do yield MathNet.Numerics.Statistics.Correlation.PearsonMatrix(path)]
                |>Seq.toArray
                |> Array.reduce(fun a b -> a+b))*1.0/float nsim
    printfn "%A" crho
    Assert.AreEqual(rho.[0,1],crho.[0,1],0.1)


[<Test>]
let TestMCMultiAssetCorrelation()=
    let nsim=1
    let ntime=1000
    let T=5.0

    let hestonpar = {   HestonParams.dividends=0.02;
                        rate=0.03;
                        kappa=16.2;
                        theta=0.01;
                        sigma=0.01;
                        rho= -0.7;
                        lambda=0.0;
                        S0=100.0;
                        V0=0.01;
                        dt=T/float ntime;
                        gamma_1=0.5;
                        gamma_2=0.5;
                        threshold_phi=1.5
                     }

    let hestonparams=Map([Heston.AssetName("A"),hestonpar;
                          Heston.AssetName("B"),hestonpar;
                          Heston.AssetName("C"),hestonpar])


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

    let spath = mcpaths.[0]
    let ld = spath
            |> Array.map(fun path ->[|1 .. path.Length-1|]
                                    |> Array.map(fun n -> path.[n]/path.[n-1]-1.0))
                            

    let rhoS = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(ld)
    printfn "Correlation rhoS \n %A" rhoS

    let rhosCube = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(cube.[0])
    printfn "Correlation rhoCube \n %A" rhosCube

    let maxdiff= (rhoS-rhosCube).ToColumnMajorArray() |> Array.maxBy(Math.Abs)
    printfn "Difference %f" maxdiff
    Assert.LessOrEqual(maxdiff,0.00001)