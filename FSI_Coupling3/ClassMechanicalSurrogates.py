from pathlib import Path

from TGaussianProcess_matern12_J7 import krig_C_7, krig_S_7, krig_W_7
from TGaussianProcess_matern12_J41 import krig_C_41, krig_S_41, krig_W_41
from TGaussianProcess_matern12_V import krig_V_C, krig_V_S, krig_V_W

BASE_DIR = Path(__file__).resolve().parent


class MechanicalSurrogate:
    """Wrapper around the three mechanical surrogate families used in the cycle."""

    def __init__(self):
        self.kriging_modelC7 = krig_C_7.KrigingModelC7(
            self._asset("TGaussianProcess_matern12_J7", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xCoefficientsC.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xiLC.npy"),
        )
        self.kriging_modelS7 = krig_S_7.KrigingModelS7(
            self._asset("TGaussianProcess_matern12_J7", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xCoefficientsS.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xiLS.npy"),
        )
        self.kriging_modelW7 = krig_W_7.KrigingModelW7(
            self._asset("TGaussianProcess_matern12_J7", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xCoefficientsW.npy"),
            self._asset("TGaussianProcess_matern12_J7", "xiLW.npy"),
        )

        self.kriging_modelVC = krig_V_C.KrigingModelVC(
            self._asset("TGaussianProcess_matern12_V", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_V", "xCoefficientsC.npy"),
            self._asset("TGaussianProcess_matern12_V", "xiLC.npy"),
        )
        self.kriging_modelVS = krig_V_S.KrigingModelVS(
            self._asset("TGaussianProcess_matern12_V", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_V", "xCoefficientsS.npy"),
            self._asset("TGaussianProcess_matern12_V", "xiLS.npy"),
        )
        self.kriging_modelVW = krig_V_W.KrigingModelVW(
            self._asset("TGaussianProcess_matern12_V", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_V", "xCoefficientsW.npy"),
            self._asset("TGaussianProcess_matern12_V", "xiLW.npy"),
        )

        self.kriging_modelC41 = krig_C_41.KrigingModelC41(
            self._asset("TGaussianProcess_matern12_J41", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xCoefficientsC.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xiLC.npy"),
        )
        self.kriging_modelS41 = krig_S_41.KrigingModelS41(
            self._asset("TGaussianProcess_matern12_J41", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xCoefficientsS.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xiLS.npy"),
        )
        self.kriging_modelW41 = krig_W_41.KrigingModelW41(
            self._asset("TGaussianProcess_matern12_J41", "xReferences.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xCoefficientsW.npy"),
            self._asset("TGaussianProcess_matern12_J41", "xiLW.npy"),
        )

    def _asset(self, *relative_parts):
        return str(BASE_DIR.joinpath(*relative_parts))

    def predict_day7(self, x_new):
        return (
            self.kriging_modelC7.kriging_function(x_new),
            self.kriging_modelS7.kriging_function(x_new),
            self.kriging_modelW7.kriging_function(x_new),
        )

    def predict_creep(self, x_new):
        return (
            30 * self.kriging_modelVC.kriging_function(x_new),
            30 * self.kriging_modelVS.kriging_function(x_new),
            30 * self.kriging_modelVW.kriging_function(x_new),
        )

    def predict_day41(self, x_new):
        return (
            self.kriging_modelC41.kriging_function(x_new),
            self.kriging_modelS41.kriging_function(x_new),
            self.kriging_modelW41.kriging_function(x_new),
        )

    # Backward-compatible wrappers kept for the existing research scripts.
    def callSurrogateModelsday7(self, xnew):
        return self.predict_day7(xnew)

    def callSurrogateModelCreep(self, xnew):
        return self.predict_creep(xnew)

    def callSurrogateModelday41(self, xnew):
        return self.predict_day41(xnew)
