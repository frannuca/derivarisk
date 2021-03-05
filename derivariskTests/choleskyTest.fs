module derivariskTests
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

    let mcpaths,cube = duration (fun () -> GBM.computeMCPaths(rho,nsim,ntime,{GBM.GBMParams.rate=0.02;dt=T/float(ntime);sigma=0.15;S0=100.0}))
    let spath = mcpaths.[0]
    let ld = spath |> Array.map(fun path ->[|1 .. path.Length-1|] |> Array.map(fun n -> path.[n]/(path.[n-1])))
    let rhoS = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(ld)    
    let rhosCube = MathNet.Numerics.Statistics.Correlation.PearsonMatrix(cube.[0])
    let drho = (rhoS-rhosCube).AsColumnMajorArray() |> Array.maxBy(Math.Abs)
    Assert.LessOrEqual(drho,0.1)

[<Test>]
let TestPriceCall()=      
    //let rho = Matrix<float>.Build.DenseOfRowArrays([|1.0|])
    let nsim=5000
    let ntime=100
    let T=0.25
    let K=90.0
    let hestonParams = {  HestonParams.dividends=0.02;
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
                          threshold_phi=1.5}
    
    let hestonparams=Map([Heston.AssetName("A"), hestonParams])

    let heston_call = Heston.heston_analytical( hestonParams.S0,
                                                K,
                                                hestonParams.kappa,
                                                0.0,
                                                T,
                                                hestonParams.rate,
                                                hestonParams.dividends,
                                                hestonParams.rho,
                                                hestonParams.sigma,
                                                hestonParams.theta,
                                                hestonParams.V0,
                                                true)

    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let assets = [Heston.AssetName "A"]
    let rho = Matrix<float>.Build.DenseOfRowArrays([|[|1.0|] |]) |> Utils.matrix2Frame assets assets

    let mcpaths, _ = duration (fun () -> Heston.computeMCPaths(rho,rho,nsim,ntime,hestonparams))
    //Utils.writeCsv(mcpaths,"/Users/fran/data/mcpath_single_asset.csv")
    
    let price= (mcpaths |> Array.map(fun assets -> assets |> Array.map(fun S -> S.[S.Length-1])|>Array.sum)
                        |> Array.map(fun ST -> Math.Max(ST-K,0.0))
                        |> Array.average)*Math.Exp(-hestonparams.[Heston.AssetName("A")].rate*T)


    printfn "Price Call Option=%f MC price =%f " heston_call price
    Assert.AreEqual(heston_call,price,0.1)

[<Test>]
let TestPricePut()=      
    //let rho = Matrix<float>.Build.DenseOfRowArrays([|1.0|])
    let nsim=5000
    let ntime=100
    let T=0.25
    let K=110.0
    let hestonParams = {  HestonParams.dividends=0.02;
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
                          threshold_phi=1.5}
    
    let hestonparams=Map([Heston.AssetName("A"), hestonParams])

    let heston_put = Heston.heston_analytical( hestonParams.S0,
                                                K,
                                                hestonParams.kappa,
                                                0.0,
                                                T,
                                                hestonParams.rate,
                                                hestonParams.dividends,
                                                hestonParams.rho,
                                                hestonParams.sigma,
                                                hestonParams.theta,
                                                hestonParams.V0,
                                                false)

    let duration f = 
        let timer = new System.Diagnostics.Stopwatch()
        timer.Start()
        let returnValue = f()
        printfn "Elapsed Time: %i" timer.ElapsedMilliseconds
        returnValue

    let assets = [Heston.AssetName "A"]
    let rho = Matrix<float>.Build.DenseOfRowArrays([|[|1.0|] |]) |> Utils.matrix2Frame assets assets

    let mcpaths, _ = duration (fun () -> Heston.computeMCPaths(rho,rho,nsim,ntime,hestonparams))
    //Utils.writeCsv(mcpaths,"/Users/fran/data/mcpath_single_asset.csv")
    
    let price= (mcpaths |> Array.map(fun assets -> assets |> Array.map(fun S -> S.[S.Length-1])|>Array.sum)
                        |> Array.map(fun ST -> Math.Max(K-ST,0.0))
                        |> Array.average)*Math.Exp(-hestonparams.[Heston.AssetName("A")].rate*T)


    printfn "Price Put Option=%f MC price =%f " heston_put price
    Assert.AreEqual(heston_put,price,0.1)


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