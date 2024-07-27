__author__ = "Rui Ribeiro"
__email__ = "rui.ribeiro@univr.it"
__year__= "2024"

#LIBRARIES
import numpy as np
import ssbtoolkit.Utils as utils

"""Module to simulate ligand-target binding curves."""
def __init__(self):
    self.receptor_conc = None
    self.lig_conc_range =  None
    self.pKd = None
    self.submax_concentration = None

def Bind(self, **kwargs):
    """
    Applies an function to calculate the fraction of occupited receptors at equilibrium.

    :parameter receptor_conc: Required (kwarg flt): concentration of receptor
    :parameter lig_conc_range: Required (kwarg array): array of range of ligand concentration
    :parameter pKd: Required (kwarg flt): pKd value of the ligand
    """

    if 'receptor_conc' not in kwargs: raise TypeError("ERROR: receptor_conc is missing")
    if 'lig_conc_range' not in kwargs: raise TypeError("ERROR: lig_conc_range is missing")
    if 'pKd' not in kwargs: raise TypeError("ERROR: pKd is missing")
    
    self._receptor_conc = kwargs.pop('receptor_conc')
    self._lig_conc_range = kwargs.pop('lig_conc_range')
    self._pKd = kwargs.pop('pKd')

    binding_data=[]
    for conc in self._lig_conc_range:
        binding_data.append(utils.CalcOccupancy(self._receptor_conc, conc, 0, self._pKd, 0))
    self.binding_data=binding_data
    return self.binding_data

def SubMaxConcentration(self):
    """
    Calculates the maximum bending point of a sigmoid-shaped curve according to the mathod of Sebaugh et al., 2003.
    
    :parameter drug_receptor: Required (int): concentration of the receptor
    :parameter lig_conc_range: Required (array): array of a range of ligand concentration
    :return: instance .submax_concentration (flt)
    
    .. note:: The minimization uses the Nelder-Mead method.
    """

    from scipy.optimize import curve_fit, minimize

    
    def sigmoid(X, Bottom, Top, Kd, p):
        return Bottom + (Top-Bottom)/(1+np.power((Kd/X),p))

    xfit = np.geomspace(np.min(self._lig_conc_range), np.max(self._lig_conc_range), 50000) #warning: shoud this be the minimum and maximum of concentration
    popt, pcov = curve_fit(sigmoid, self._lig_conc_range, self.binding_data, bounds=([np.min(self.binding_data),-np.inf,-np.inf, 0.5],[np.inf,np.max(self.binding_data),np.inf, 2.5]))


    def sigmoid_deriv_b(x, a,d,c,b):
        return (x/c)**b*(a - d)*np.log(x/c)/((x/c)**b + 1)**2

    min_value = minimize(sigmoid_deriv_b, np.max(xfit), args=(popt[0],popt[1],popt[2],popt[3]), method = 'Nelder-Mead')

    self.submax_concentration = round(min_value.x[0],3)
    return self.submax_concentration

def ShowCurve(self):
    """
    Plots ligand-target binding curve
    """

    #import plotly
    import plotly.graph_objs as go 
    from scipy.optimize import curve_fit
    from sklearn.preprocessing import minmax_scale

    ##Fitting curve to the data

    yy = minmax_scale(self.binding_data)*100

    def equation_dose(X, Bottom, Top, EC50, p):
        return Bottom + (Top-Bottom)/(1+np.power((EC50/X),p))

    popt, pcov = curve_fit(equation_dose, self._lig_conc_range, yy, bounds=([np.min(yy),-np.inf,-np.inf, 0.5],[np.inf,np.max(yy),np.inf, 2.5]))

    xfit = np.geomspace(np.min(self._lig_conc_range), np.max(self._lig_conc_range), 50000) # These values are the same as the values for the simulation time and not ligand concentration
    yfit = minmax_scale(equation_dose(xfit, *popt))*100

    

    trace1 = go.Line(x=xfit, y=yfit, showlegend=False, name='radioligand')
    if self.submax_concentration:
        xsubmaximal = np.array(self.submax_concentration)
        ysubmaximal = np.array(equation_dose(xsubmaximal, *popt))
        trace2 = go.Scatter(x=xsubmaximal, y=ysubmaximal, showlegend=True, mode='markers', name='submaximal ({} μM)'.format(xsubmaximal),
                    marker=dict(size=14))
    else: trace2=[]


    layout = dict(title = '',
                            xaxis = dict(
                                title = '[ligand] μM',
                                type ='log',
                                exponentformat='e',
                                titlefont=dict(
                                    size=20
                                ),
                                tickfont=dict(
                                    size=20
                                )),
                            yaxis = dict(
                                title = '% occupied receptors',
        
                                titlefont=dict(
                                    size=20),
                                tickfont=dict(
                                size=20)

                            ),
                            legend=dict(font=dict(size=15)),
                            autosize=False,
                            width=850,
                            height=650,
                            )
    if self.submax_concentration:
        fig = go.Figure(data=[trace1, trace2], layout=layout)
    else:
        fig = go.Figure(data=[trace1], layout=layout)
    return fig
