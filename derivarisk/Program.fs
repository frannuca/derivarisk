// Learn more about F# at http://fsharp.org

open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
[<EntryPoint>]
let main argv =

    
    let rs=0.95
    //let rho = Matrix<float>.Build.DenseOfRowArrays(
    //    [|
    //        [|1.0; rs;rs|];
    //        [|rs; 1.0;rs|];
    //        [|rs; rs;1.0|]
    //    |]
    //)
    let rho = Matrix<float>.Build.DenseOfRowArrays([|1.0|])
    let nsim=5000
    let ntime=100
    let T=0.25
    //let hestonparams=Map([Heston.AssetName("A"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.005;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5};
    //                      Heston.AssetName("B"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.005;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5};
    //                      Heston.AssetName("C"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.005;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5}
    //                     ])

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
    let mrho = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(mcpaths.[10])
    printfn "Correlation %A" mrho
    
    0 // return an integer exit code
