import numpy as np
from TGaussianProcess_matern12_J7 import krig_C_7, krig_S_7, krig_W_7 
from TGaussianProcess_matern12_V import krig_V_C, krig_V_S, krig_V_W
from TGaussianProcess_matern12_J41 import krig_C_41, krig_S_41, krig_W_41

class MechanicalSurrogate:
    
    
    def __init__(self):
        
        self.kriging_modelC7 = krig_C_7.KrigingModelC7('TGaussianProcess_matern12_J7/xReferences.npy', 'TGaussianProcess_matern12_J7/xCoefficientsC.npy', 
                                          'TGaussianProcess_matern12_J7/xiLC.npy')
        self.kriging_modelS7 = krig_S_7.KrigingModelS7('TGaussianProcess_matern12_J7/xReferences.npy', 'TGaussianProcess_matern12_J7/xCoefficientsS.npy', 
                                          'TGaussianProcess_matern12_J7/xiLS.npy')
        self.kriging_modelW7 = krig_W_7.KrigingModelW7('TGaussianProcess_matern12_J7/xReferences.npy', 'TGaussianProcess_matern12_J7/xCoefficientsW.npy', 
                                          'TGaussianProcess_matern12_J7/xiLW.npy')
        self.kriging_modelVC = krig_V_C.KrigingModelVC('TGaussianProcess_matern12_V/xReferences.npy', 'TGaussianProcess_matern12_V/xCoefficientsC.npy', 
                                          'TGaussianProcess_matern12_V/xiLC.npy')
        self.kriging_modelVS = krig_V_S.KrigingModelVS('TGaussianProcess_matern12_V/xReferences.npy', 'TGaussianProcess_matern12_V/xCoefficientsS.npy', 
                                          'TGaussianProcess_matern12_V/xiLS.npy')
        self.kriging_modelVW = krig_V_W.KrigingModelVW('TGaussianProcess_matern12_V/xReferences.npy', 'TGaussianProcess_matern12_V/xCoefficientsW.npy', 
                                          'TGaussianProcess_matern12_V/xiLW.npy')
        
        self.kriging_modelC41 = krig_C_41.KrigingModelC41('TGaussianProcess_matern12_J41/xReferences.npy', 'TGaussianProcess_matern12_J41/xCoefficientsC.npy', 
                                          'TGaussianProcess_matern12_J41/xiLC.npy')
        self.kriging_modelS41 = krig_S_41.KrigingModelS41('TGaussianProcess_matern12_J41/xReferences.npy', 'TGaussianProcess_matern12_J41/xCoefficientsS.npy', 
                                          'TGaussianProcess_matern12_J41/xiLS.npy')
        self.kriging_modelW41 = krig_W_41.KrigingModelW41('TGaussianProcess_matern12_J41/xReferences.npy', 'TGaussianProcess_matern12_J41/xCoefficientsW.npy', 
                                          'TGaussianProcess_matern12_J41/xiLW.npy')
    
    def callSurrogateModelsday7(self,xnew):
        
        c7 = self.kriging_modelC7.kriging_function(xnew)
        s7 = self.kriging_modelS7.kriging_function(xnew)
        w7 = self.kriging_modelW7.kriging_function(xnew)
        
        return c7,s7,w7
    
    def callSurrogateModelCreep(self,xnew):

        c_creep = 30 * self.kriging_modelVC.kriging_function(xnew)
        s_creep = 30 * self.kriging_modelVS.kriging_function(xnew)
        w_creep = 30 * self.kriging_modelVW.kriging_function(xnew)

        return c_creep , s_creep, w_creep
    
    def callSurrogateModelday41(self,xnew):

        c41 = self.kriging_modelC41.kriging_function(xnew)
        s41 = self.kriging_modelS41.kriging_function(xnew)
        w41 = self.kriging_modelW41.kriging_function(xnew)
        
        return c41,s41,w41



































"""x_0 = np.array([-0.00113202 ,-0.00105823 ,0.00016258 ,-115 ,-50.0, -30 ,-40.00000000 ,549.5870518112287, 8.361867242525591e+17 ,0.01993624709493367, 1, 0.8117634667679304 ,1.3738148268682195 ])
s = MechanicalSurrogate().callSurrogateModelday41(x_0)

print(s)"""

