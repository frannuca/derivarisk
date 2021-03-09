namespace DerivaRisk
open System
open MathNet.Numerics.LinearAlgebra

type OptionType=
    |CALL
    |PUT


type BARRIERTYPE=
    |DownOut of float
    |UpOut of float
    |DownIn of float
    |UpIn of float

[<Measure>] type years
[<Measure>] type months
[<Measure>] type days


[<Interface>]
type IOptionDataBase=
     abstract StartDate:DateTime  with get
     abstract EndDate:DateTime with get
     abstract SettlementDate:DateTime with get     
     
type Spot=Spot of float

[<Interface>]
type MarketDataBase=
    abstract rate:(float<years> -> float) with get
    abstract DividendYield:DateTime option with get
    abstract ProportionalYield:Map<DateTime,float> option with get
    abstract volatility: Spot -> float<years> -> float

type SimCubeData={nsim:int;ntimesteps:int;seed:int;number_assets:int;rho:Matrix<float>}
    