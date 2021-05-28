###################################
#To regenerate instrument classes:#
###################################

xsd *.xsd /c /o:.. /n:InstrumentData
mv ../commons_instrument_marketdata_portfolio_simulation.cs ../messagedata.cs