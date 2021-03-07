// Learn more about F# at http://fsharp.org

open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Heston
open GBM
open Deedle

[<EntryPoint>]
let main argv =
    let nsim=20000
    let ntime=50
    let T=0.25
    let hestonParams = {    HestonParams.dividends=0.02;
                            rate=0.03;
                            kappa=6.2;
                            theta=0.06;
                            sigma=0.20;
                            rho= -0.7;
                            lambda=0.0;
                            S0=100.0;
                            V0=0.03;
                            dt=T/float ntime;
                            gamma_1=0.5;
                            gamma_2=0.5;
                            threshold_phi=1.5}

    
   
    let rhos = [|[|1.0|]|]
                |> Matrix<float>.Build.DenseOfColumnArrays
                |> Utils.matrix2Frame [Heston.AssetName "A"] [Heston.AssetName "A"]

    let getpar x=
        [(Heston.AssetName("A"),{hestonParams with S0=x})]|>Map.ofSeq
    let Smin=90.0
    let Smax=125.0
    let nS = 200
    let K=100.0
    let barrier = HestonBarrier.BARRIERTYPE.UpOut(120.0)

    //let aaa= HestonBarrier.barrierprice(nsim,ntime,getpar Smin,rhos,rhos,K,HestonBarrier.BARRIERTYPE.DownOut(90.0),Utils.CALL)
    //let prices_callupout = Utils.linspace Smin Smax nS
    //                        |> Array.map(fun S->
    //                                printfn "Simulating S=%f in range [%f,%f]" S Smin Smax
    //                                S,
    //                                HestonBarrier.barrierprice(nsim,ntime,getpar S,rhos,rhos,K,HestonBarrier.BARRIERTYPE.UpOut(110.0),Utils.CALL))                    
    //                        |> Series.ofObservations
    //let prices_calldownout = Utils.linspace Smin Smax nS
    //                            |> Array.map(fun S->
    //                                    printfn "Simulating S=%f in range [%f,%f]" S Smin Smax
    //                                    S,
    //                                    HestonBarrier.barrierprice(nsim,ntime,getpar S,rhos,rhos,K,HestonBarrier.BARRIERTYPE.DownOut(85.0),Utils.CALL))                    
    //                            |> Series.ofObservations



    //let prices_callupin = Utils.linspace 30.0 80.0 nS
    //                        |> Array.map(fun S->
    //                                printfn "Simulating S=%f in range [%f,%f]" S Smin Smax
    //                                S,
    //                                HestonBarrier.barrierprice(nsim,ntime,getpar S,rhos,rhos,40.0,HestonBarrier.BARRIERTYPE.UpIn(65.0),Utils.CALL))                    
    //                        |> Series.ofObservations

    let saxis = Utils.linspace Smin Smax nS
    let ds = saxis.[1]-saxis.[0]
    
    let mcprice = saxis
                            |> Array.Parallel.map(fun S->
                                    printfn "Simulating S=%f in range [%f,%f]" S Smin Smax
                                    S,
                                    HestonBarrier.barrierprice(nsim,ntime,getpar S,rhos,rhos,K,barrier,Utils.CALL))                    
                            |> Series.ofObservations |>Series.sortByKey

    let prices_call = saxis
                              |> Array.Parallel.map(fun S->
                                      printfn "Simulating S=%f in range [%f,%f]" S Smin Smax
                                      S,
                                      Heston.heston_analytical(S,100.0,hestonParams.kappa,0.0,T,hestonParams.rate,hestonParams.dividends,hestonParams.rho,hestonParams.sigma,hestonParams.theta,hestonParams.V0,true)                    )
                              |> Series.ofObservations |>Series.sortByKey                                
    printfn("FINISHED  Monte Carlo Calculation")
    let deltas = (mcprice |> Series.diff(1)) / ds
    let gammas =  (deltas |> Series.diff(1)) / ds
    printfn("FINISHED  Greeks Calculation")
    //let frame = Frame.ofColumns["UpOut"=>prices_callupout;"Call_DownOut"=>prices_calldownout;"Call_UpIn"=>prices_callupin;"Call_DownIn"=>prices_calldownIn]
    let frame = Frame.ofColumns["analytical"=>prices_call;"prices_call_downIn"=>mcprice;"delta"=>deltas;"gamma"=>gammas]
    printfn("Creating Frame")
    frame.SaveCsv("/Users/fran/data/price_analysis.csv",["Spot"],',')
    printfn("Frame saved")
    0 // return an integer exit code
