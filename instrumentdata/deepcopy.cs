using System;
using System.Xml;
using System.Xml.Serialization;
using System.IO;

namespace InstrumentData
{ 

    static public class deepcopy
    {
        public static T DeepCopy<T>(T other)
        {
            var xmlSerializer = new XmlSerializer(typeof(T));

            string xmlstr = "";
            using (StringWriter textWriter = new StringWriter())
            {
                xmlSerializer.Serialize(textWriter, other);
                xmlstr= textWriter.ToString();                
            }
            using(StringReader textReader = new StringReader(xmlstr))
            {
                return (T)xmlSerializer.Deserialize(textReader);
            }
        }

        static public Marketdata shockMarketData(Marketdata md,double shock)
        {
            var md_p = DeepCopy(md);
            foreach(var s in md_p.spots)
            {
                foreach(var spot in s.points)
                {
                    spot.y *= (1.0 + shock);
                }
            }

            return md_p;
        }
    }

    
}
