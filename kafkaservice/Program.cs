// Copyright 2020 Confluent Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Refer to LICENSE for more information.

using Confluent.Kafka;
using Confluent.Kafka.SyncOverAsync;
using System;
using System.ComponentModel.DataAnnotations;
using System.Threading;
using System.Threading.Tasks;
using Newtonsoft.Json;
using InstrumentData;

/// <summary>
///     An example of working with JSON data, Apache Kafka and 
///     Confluent Schema Registry (v5.5 or later required for
///     JSON schema support).
/// </summary>
namespace Confluent.Kafka.Examples.JsonSerialization
{
    /// <summary>
    ///     A POCO class corresponding to the JSON data written
    ///     to Kafka, where the schema is implicitly defined through 
    ///     the class properties and their attributes.
    /// </summary>
    /// <remarks>
    ///     Internally, the JSON serializer uses Newtonsoft.Json for
    ///     serialization and NJsonSchema for schema creation and
    ///     validation. You can use any property annotations recognised
    ///     by these libraries.
    ///
    ///     Note: Off-the-shelf libraries do not yet exist to enable
    ///     integration of System.Text.Json and JSON Schema, so this
    ///     is not yet supported by the Confluent serializers.
    /// </remarks>
    class Person
    {
        [JsonRequired] // use Newtonsoft.Json annotations
        [JsonProperty("firstName")]
        public string FirstName { get; set; }

        [JsonRequired]
        [JsonProperty("lastName")]
        public string LastName { get; set; }

        [Range(0, 150)] // or System.ComponentModel.DataAnnotations annotations
        [JsonProperty("age")]
        public int Age { get; set; }
    }

    class Program
    {
        static async Task Main(string[] args)
        {
            //if (args.Length != 3)
            //{
            //    Console.WriteLine("Usage: .. bootstrapServers schemaRegistryUrl topicName");
            //    return;
            //}

            string bootstrapServers = "localhost:9092";// args[0];          
            string topicName = "analytics"; //args[2];

            var producerConfig = new ProducerConfig
            {
                BootstrapServers = bootstrapServers
            };

           
            var consumerConfig = new ConsumerConfig
            {
                BootstrapServers = bootstrapServers,
                GroupId = "json-example-consumer-group"
            };
          

            CancellationTokenSource cts = new CancellationTokenSource();
            var consumeTask = Task.Run(() =>
            {
                using (var consumer =
                    new ConsumerBuilder<Null, string>(consumerConfig)                       
                        .SetErrorHandler((_, e) => Console.WriteLine($"Error: {e.Reason}"))
                        .Build())
                {
                    consumer.Subscribe(topicName);

                    try
                    {
                        while (true)
                        {
                            try
                            {
                                var cr = consumer.Consume(cts.Token);

                                var amount = InstrumentSerializers.String2Instrument<InstrumentData.Amount>(cr.Message.Value);

                                //Console.WriteLine($"Name: {person.FirstName} {person.LastName}, age: {person.Age}");
                                Console.WriteLine($"Amount: {amount.amount} {amount.amountType}, type: {amount.amountType}");
                            }
                            catch (ConsumeException e)
                            {
                                Console.WriteLine($"Consume error: {e.Error.Reason}");
                            }
                        }
                    }
                    catch (OperationCanceledException)
                    {
                        consumer.Close();
                    }
                }
            });


            using (var producer = new derivarisk.fkafkaservice.Producer<Amount>(new derivarisk.fkafkaservice.ProducerConfig("localhost:9092", "analytics")))
            {
                for (int i = 0; i < 100; ++i)
                {
                    var r =
                    producer.send(new Amount()
                    {
                        amount = i,
                        amountType = AmountType.PERCENTAGE
                    });
                    
                    Console.WriteLine($"Send amount {i} with result {r.Status}");
                }
            }
               
            

            cts.Cancel();

            
        }
    }
}