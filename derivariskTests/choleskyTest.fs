module derivariskTests
open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open NUnit.Framework

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
let TestPriceCall()=
    
    let rho = Matrix<float>.Build.DenseOfRowArrays(
        [|
            [|1.0|]
        |]
    )
    //let rho = Matrix<float>.Build.DenseOfRowArrays([|1.0|])
    let nsim=5000
    let ntime=100
    let T=0.25
    
    let hestonparams=Map([Heston.AssetName("A"),
        {HestonParams.dividends=0.02;
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
        }])

    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let mcpaths = duration (fun () -> Heston.computeMCPaths(rho,nsim,ntime,hestonparams))
    let K=90.0
    let price= (mcpaths |> Array.map(fun assets -> assets |> Array.map(fun S -> S.[S.Length-1])|>Array.sum)
                        |> Array.map(fun ST -> Math.Max(ST-K,0.0))
                        |> Array.average)*Math.Exp(-hestonparams.[Heston.AssetName("A")].rate*T)

    printfn "Price Call Option=%f" price
    Assert.AreEqual(11.26,price,0.5)