module HestonMultiAssetTest
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
let TestMCMultiAssetPrice()=
    let rs=0.95
    let assets = [Heston.AssetName "A";Heston.AssetName "B";Heston.AssetName "C"]
    let rho = Matrix<float>.Build.DenseOfRowArrays(
                            [|  [|1.0; rs;rs|];
                                [|rs; 1.0;rs|];
                                [|rs; rs;1.0|]
                            |]) |> Utils.matrix2Frame assets assets
    let rv = rs
    let rho_vol = Matrix<float>.Build.DenseOfRowArrays(
                        [|
                            [|1.0; rv;rv|];
                            [|rv; 1.0;rv|];
                            [|rv; rv;1.0|]
                        |])|> Utils.matrix2Frame assets assets

    let nsim=50
    let ntime=10
    let T=0.25
    let hestonparams=Map([Heston.AssetName("A"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.5;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5};
                          Heston.AssetName("B"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.5;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5};
                          Heston.AssetName("C"), {HestonParams.dividends=0.02;rate=0.03;kappa=6.2;theta=0.6;sigma=0.5;rho= -0.7;lambda=0.0;S0=100.0;V0=0.03;dt=1.0/float ntime;gamma_1=0.5;gamma_2=0.5;threshold_phi=1.5}
                        ])

    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let mcpaths, _ = duration (fun () -> Heston.computeMCPaths(rho,rho_vol,nsim,ntime,hestonparams))
    use stream = new StreamWriter("/Users/fran/data/mcpath.csv",false)
    mcpaths |> Array.iteri(fun nsim assets -> assets|>Array.iter(fun asset -> stream.WriteLine(nsim.ToString()+","+String.Join(",",asset))))

    let K=90.0
    let price= (mcpaths |> Array.map(fun assets -> assets |> Array.map(fun S -> S.[S.Length-1])|>Array.sum)
                        |> Array.map(fun ST -> Math.Max(ST-K,0.0))
                        |> Array.average)*Math.Exp(-hestonparams.[Heston.AssetName("A")].rate*T)

    printfn "Price Call Option=%f" price
    let mrho = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(mcpaths.[10])
    printfn "Correlation %A" mrho