namespace DerivaRisk.model
open InstrumentData
open DerivaRisk.data

[<AbstractClass>]
type MonteCarloModel(simconfig:SimulationCube,instr:Instrument,md:Marketdata)=
    let _config=simconfig
    let _instr=instr
    let _md=md
    

    abstract generatepath: unit->float[][][]
    abstract payoff:float[][]->float
    
    interface IModel with
        member self.compute(flags:ComputationSelection array)
                           (instr:Instrument)
                           (md:Marketdata):ModelResult=


            let cube = self.generatepath()
            let npv = cube |> Array.map(self.payoff) |> Array.average

            let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            let T = float(instr.expiry.Subtract(md.referenceDate).TotalDays)/360.0

            let r = Extractors.Integrate_Rate md instr.currency (t,T)
            let q = Extractors.Integrate_Dividends md instr.currency (t,T)
            

            {ModelResult.npv=npv*exp(-(r-q)*(T-t));delta=None;gamma=None;theta=None;rho=None}

    member val instr = _instr
    member val md = _md
    member val simconfig = _config
    member val optionbase = instr :?> OptionBase