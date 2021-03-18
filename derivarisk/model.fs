namespace DerivaRisk.model

open InstrumentData
open MathNet.Numerics.LinearAlgebra

type ModelResult={npv:float;delta:Vector<float> option;gamma:Matrix<float> option;theta:float option;rho:float option}
type ComputationSelection=
    |NPV
    |Delta
    |Gamma
    |Theta
    |Rho

[<Interface>]
type IModel=
    abstract compute: ComputationSelection array->Instrument->Marketdata->ModelResult
