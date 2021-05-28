using System;
using Newtonsoft.Json;

namespace kafkaservice
{
    static public class xmlMessageConverter
    {
        public static string Xml2String<T>(T obj)
        {
            return JsonConvert.SerializeObject(obj);
        }

        public static T String2XmlObj<T>(string xml)
        {
            T obj = JsonConvert.DeserializeObject<T>(xml);
            return obj;
        }


    }
}
