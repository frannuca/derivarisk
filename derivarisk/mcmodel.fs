namespace DerivaRisk.model
open InstrumentData
open DerivaRisk.data

[<AbstractClass>]
type MonteCarloModel(simconfig:SimulationCube,instr:Instrument)=
    let _config=simconfig
    let _instr=instr
    

    abstract generatepath: Marketdata->float[][][]
    abstract payoff:float[][]->float
    
    interface IModel with
        member self.npv(instr:Instrument)
                       (md:Marketdata):float=


            let cube = self.generatepath md
            let npv = cube |> Array.map(self.payoff) |> Array.average

            let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            let T = float(instr.expiry.Subtract(md.referenceDate).TotalDays)/360.0

            let r = Extractors.Integrate_Rate md instr.currency (t,T)
            let q = Extractors.Integrate_Dividends md instr.currency (t,T)            
            npv

        member self.greekDelta(instr:Instrument)
                              (md:Marketdata):float array =
                      
            let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            let T = float(instr.expiry.Subtract(md.referenceDate).TotalDays)/360.0

            let vanillaoption = instr :?> OptionBase
            let md_p = deepcopy.shockMarketData(md,0.01)

            let model = self :> IModel
            let npv = model.npv instr md

            let compute_internal_delta underlyingId =                               
                let npv_p = model.npv instr md_p
                let So = Extractors.extract_spot_for_time_t md t underlyingId 
                let Sp = Extractors.extract_spot_for_time_t md_p t underlyingId
                (npv_p-npv)/(Sp-So)

            vanillaoption.basket.components
            |> Array.map(compute_internal_delta)
            
        member self.greekGamma(instr:Instrument)
                              (md:Marketdata):float array =
              
            let t = float(instr.startdate.Subtract(md.referenceDate).TotalDays)/360.0
            
            let vanillaoption = instr :?> OptionBase
            let h=0.01
            let md_p = deepcopy.shockMarketData(md,h)
            let md_m = deepcopy.shockMarketData(md,-h)

            let model = self :> IModel
            let npv = model.npv instr md
            let npv_p = model.npv instr md_p
            let npv_m = model.npv instr md_m

            let compute_internal_delta underlyingId =                                               
                let So = Extractors.extract_spot_for_time_t md t underlyingId 
                let Sp = Extractors.extract_spot_for_time_t md_p t underlyingId               
                (npv_p-2.0*npv+npv_m)/(Sp-So)**2.0

            vanillaoption.basket.components
            |> Array.map(compute_internal_delta)


    member val instr = _instr    
    member val simconfig = _config
    member val optionbase = instr :?> OptionBase