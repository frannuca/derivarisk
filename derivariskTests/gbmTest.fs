module GBMTests
open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open NUnit.Framework
open GBM
open System.IO

[<SetUp>]
let Setup () =
    ()

[<Test>]
let TestGBMCorrelation()=
    let rs= 0.6
    let rho = Matrix<float>.Build.DenseOfRowArrays(
                                        [|
                                            [|1.0; rs;rs|];
                                            [|rs; 1.0;rs|];
                                            [|rs; rs;1.0|]
                                        |])

    let nsim=1
    let ntime=1000
    let T=1.0

    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let par = {GBMParams.rate=0.02;dt=T/float(ntime);sigma= (fun _ _ -> 0.15);S0=100.0}

    let mcpar = {SimCubeData.nsim=nsim;
                 SimCubeData.ntimesteps=ntime;
                 SimCubeData.number_assets=rho.ColumnCount;
                 SimCubeData.rho=rho;
                 SimCubeData.seed=seedS;
                }

    let mcpaths,cube = duration (fun () -> GBM.computeMCPaths mcpar (Array.create (rho.ColumnCount) par))
    
    let spath = mcpaths.[0]
    let ld = spath |> Array.map(fun path ->[|1 .. path.Length-1|] |> Array.map(fun n -> path.[n]/(path.[n-1])))
    let rhoS = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(ld)    
    let rhosCube = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(cube.[0])
    let drho = (rhoS-rhosCube).AsColumnMajorArray() |> Array.maxBy(Math.Abs)
    Assert.LessOrEqual(drho,0.1)
