namespace derivarisk.fkafkaservice

open Confluent.Kafka;
open Confluent.Kafka.SyncOverAsync;
open System;
open System.ComponentModel.DataAnnotations;
open System.Threading;
open System.Threading.Tasks;
open Newtonsoft.Json;
open System.IO
open System
open InstrumentData
open System.Net

type ProducerConfig={ Server:string; Topic:string}

type Producer<'T>(config:ProducerConfig)=
    
    let producer_config =
         let p = ProducerConfig()
         p.BootstrapServers <- config.Server
         p


    let _producer = ProducerBuilder<Null,string>(producer_config).Build()


    interface IDisposable with 
        member self.Dispose() =
            _producer.Dispose()


    member self.send(obj:'T)=        
            let msg= new Message<Null,string>()
            msg.Value <- InstrumentData.InstrumentSerializers.Instrument2String<'T>(obj)
            _producer.ProduceAsync(config.Topic,msg).Result
            

        
    member self.sendAsync(obj:'T)=
        async{
            
            let msg= new Message<Null,string>()
            msg.Value <- InstrumentData.InstrumentSerializers.Instrument2String<'T>(obj)            
            let! r = Async.AwaitTask(_producer.ProduceAsync(config.Topic,msg))
            return r
        }
        