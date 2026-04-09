import ROOT

# Create shortcuts if uranie exists
urasys = ROOT.TString(ROOT.gSystem.Getenv("URANIESYS"))
if not urasys.EqualTo(""):
    from ROOT.URANIE import DataServer as DataServer
    from ROOT.URANIE import Sampler as Sampler
    from ROOT.URANIE import Launcher as Launcher
    from ROOT.URANIE import Relauncher as Relauncher
    from ROOT.URANIE import Reoptimizer as Reoptimizer
    from ROOT.URANIE import Sensitivity as Sensitivity
    from ROOT.URANIE import Optimizer as Optimizer
    from ROOT.URANIE import Modeler as Modeler
    from ROOT.URANIE import Calibration as Calibration
    from ROOT.URANIE import UncertModeler as UncertModeler
    from ROOT.URANIE import Reliability as Reliability
    from ROOT.URANIE import XMLProblem as XMLProblem
    from ROOT.URANIE import MpiRelauncher as MpiRelauncher
    pass

# General graphical style
white = 0

# PlotStyle
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptDate(21)

# Legend
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetFillStyle(0)

# Pads
ROOT.gStyle.SetPadColor(white)
ROOT.gStyle.SetTitleFillColor(white)
ROOT.gStyle.SetStatColor(white)

#  ====================  Hint ====================
#
#    Might be practical to store this in a convenient place (for instance
#    the ".python" folder in your home directory) or any other place where
#    your $PYTHONPATH is pointing.
#
#    example : export PYTHONPATH=$PYTHONPATH:${HOME}/.mypython/
#
#    It should then be called as "from rootlogon import " + the list of module
#    This would replace the shortcuts created and import done in the rest of
#    the scripts
#
#    Many style issue can be set once and for all here.
#    toto=DataServer.TDataServer()
#
