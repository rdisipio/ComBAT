import PyCintex
PyCintex.loadDict("libCrossSectionFinder.so")
from ROOT import CrossSectionFinder, BCLog, BCAux, BCDataSet, BCDataPoint

BCAux.SetStyle()
BCLog.SetLogLevel(BCLog.detail);

m = CrossSectionFinder()

backgroundMeasurement = BCDataPoint(2)
backgroundMeasurement.SetValue(0, 100)
backgroundMeasurement.SetValue(1, 100)

dataSet = BCDataSet()
BCDataSet.AddDataPoint(backgroundMeasurement)

m.SetDataSet(dataSet)

BCLog.OutSummary("Test model created")

m.MarginalizeAll()
m.FindMode( m.GetBestFitParameters() )

m.PrintResults("combination_results.txt")
BCLog.CloseLog()

BCLog.OutSummary("Test program ran successfully")

