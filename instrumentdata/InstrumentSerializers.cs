using System;
using System.Xml;
using System.Xml.Serialization;
using System.IO;

namespace InstrumentData
{
    static public class InstrumentSerializers
    {
        public static string Instrument2String<T>(T obj)
        {
            var xmlSerializer = new XmlSerializer(typeof(T));

            string xmlstr = "";
            using (StringWriter textWriter = new StringWriter())
            {
                xmlSerializer.Serialize(textWriter, obj);
                xmlstr = textWriter.ToString();
            }
            return xmlstr;
        }

        public static T String2Instrument<T>(string xmlstr)
        {
            var xmlSerializer = new XmlSerializer(typeof(T));           
            using (StringReader textReader = new StringReader(xmlstr))
            {
                return (T)xmlSerializer.Deserialize(textReader);
            }
        }
    }
}
