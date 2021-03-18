// Learn more about F# at http://fsharp.org

open System
open DerivaRisk
open MathNet.Numerics.LinearAlgebra
open Deedle
open InstrumentData
open DerivaRisk.model

[<EntryPoint>]
let main argv =

    let basket = Basket()
    basket.components <- Array.init(1)(fun _ -> Identifier())

    basket.components.[0]<- Identifier()
    basket.components.[0].id<-"AAA"
    basket.components.[0].name<-"Stock"
    basket.components.[0].weight<-1.0



    let voption = OptionBase()
    voption.basket<-basket;

    voption.currency<- Currency()
    voption.currency.id<-"chf"
    voption.currency.isocode<-"CHF"
    voption.currency.name<-"Swiss Franc"

    voption.expiry <- DateTime(2020,03,31)
    voption.startdate <- DateTime(2020,01,01)

    voption.id <- "European Call"
    voption.optionType<-OptionType.CALL
    voption.strike<-Amount()
    voption.strike.amount<-90.0
    voption.strike.amountType <- AmountType.ABSOLUTE

    let simconfig = SimulationCube()
    simconfig.number_of_assets<- uint32(1)
    simconfig.number_of_simulations <- uint32(5000)
    simconfig.number_of_time_steps <- uint32(100)
    simconfig.seed <- uint32(42)
    simconfig.seed2 <- uint32(51)
    simconfig.dt <- float32(0.25/float simconfig.number_of_time_steps)
    simconfig.dt1 <- float32(1.0/365.0)
    simconfig.assets_correlation <- CorrelationMatrix()
    simconfig.assets_correlation.assets <- [|"AAA"|]
    simconfig.assets_correlation.isRowMajor<-true
    simconfig.assets_correlation.data <- [|1.0|]


    let md = Marketdata()
    md.rates <- Array.init 1 (fun _ -> MDseries())
    md.rates.[0].id<-"chf"
    md.rates.[0].points <- [|Point();Point()|]
    md.rates.[0].points.[0].x <- 0.0
    md.rates.[0].points.[0].y <- 0.03
    md.rates.[0].points.[1].x <- 1.0
    md.rates.[0].points.[1].y <- 0.03
    md.rates.[0].timeunits <- TimeUnits.years

    md.dividends <- Array.init 1 (fun _ -> MDseries())
    md.dividends.[0].id<-"chf"
    md.dividends.[0].points <- [|Point();Point()|]
    md.dividends.[0].points.[0].x <- 0.0
    md.dividends.[0].points.[0].y <- 0.02
    md.dividends.[0].points.[1].x <- 1.0
    md.dividends.[0].points.[1].y <- 0.02
    md.dividends.[0].timeunits <- TimeUnits.years

    md.stochasticVolatilities <- [|HestonVolatility()|]
    md.stochasticVolatilities.[0].gamma1 <- 0.5
    md.stochasticVolatilities.[0].gamma2 <- 0.5
    md.stochasticVolatilities.[0].id <- "AAA"
    md.stochasticVolatilities.[0].kappa <- 6.2
    md.stochasticVolatilities.[0].lambda <- 0.0
    md.stochasticVolatilities.[0].name <- "AAA"
    md.stochasticVolatilities.[0].rho <- -0.7 
    md.stochasticVolatilities.[0].sigma <- 0.5
    md.stochasticVolatilities.[0].theta <- 0.06
    md.stochasticVolatilities.[0].Vo <- 0.03

    md.referenceDate <- DateTime(2020,01,01)

    md.spots <- [|MDseries()|]
    md.spots.[0].id<-"AAA"
    md.spots.[0].weight <- 1.0
    md.spots.[0].points <- [|Point();Point()|]
    md.spots.[0].points.[0].x <- 0.0
    md.spots.[0].points.[0].y <- 100.0
    md.spots.[0].points.[1].x <- 1.0
    md.spots.[0].points.[1].y <- 100.0
    md.spots.[0].timeunits <- TimeUnits.years
    let pricer= MonteCarloHeston(simconfig,voption,md) :> IModel
    let npv = pricer.compute([|ComputationSelection.NPV|]) voption md

    let analytical = HestonAnalytical() :> IModel
    let npv2 = analytical.compute([|ComputationSelection.NPV|]) voption md
    0 // return an integer exit code
    (*
    HestonParams.dividends=0.02;
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
    threshold_phi=1.5*)