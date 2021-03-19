namespace DerivaRisk.model

open InstrumentData
open MathNet.Numerics.LinearAlgebra

type ModelResult={npv:float;delta:float array option;gamma:Matrix<float> option;theta:float option;rho:float option}
type ComputationSelection=
    |NPV
    |Delta
    |Gamma
    |Theta
    |Rho

[<Interface>]
type IModel=
    abstract npv: Instrument->Marketdata->float
    abstract greekDelta: Instrument->Marketdata->float array
    abstract greekGamma: Instrument->Marketdata->float array
