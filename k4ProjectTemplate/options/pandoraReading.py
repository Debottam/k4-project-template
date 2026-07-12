from Configurables import pandoraReading
from k4FWCore import ApplicationMgr, IOSvc
from Gaudi.Configuration import DEBUG

io = IOSvc(
    Input=["input_pandora.root"],
    Output="output_pandora.root",
    CollectionNames=["PandoraPFOs"]
)

alg = pandoraReading(
    "PandoraReading",
    InputPFOs="PandoraPFOs",
    OutputPFOs="PandoraPFOsOutput"
)

ApplicationMgr(
    TopAlg=[alg],
    ExtSvc=[io],
    EvtSel="NONE",
    EvtMax=2,
    OutputLevel=DEBUG,
)
