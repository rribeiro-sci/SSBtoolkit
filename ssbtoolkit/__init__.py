import pkgutil

__all__ = []
for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    __all__.append(module_name)
    _module = loader.find_module(module_name).load_module(module_name)
    globals()[module_name] = _module


import sys, os, warnings
sys.path.insert(0, os.path.abspath(os.path.split(os.path.realpath(__file__))[0])) 
warnings.simplefilter(action='ignore')




import platform, site 
distpath = site.getsitepackages()[0]
if platform.system() == 'Linux':
    BioNetGen=os.path.join(distpath, 'bionetgen/bng-linux:')
elif platform.system() == 'Darwin':
    BioNetGen=os.path.join(distpath, 'bionetgen/bng-mac:')
elif platform.system()=='Windows':
    BioNetGen=os.path.join(distpath, 'bionetgen/bng-win:')
else:
    raise ValueError('BioNetGen error. Platform unknown! The pygomodo was tested in Linux and Darwin (Mac) platforms.')
os.environ['PATH']=BioNetGen+os.environ['PATH']


#LIBRARIES
import importlib
import pylab as pl
from pysb.simulator import ScipyOdeSimulator
from scipy.optimize   import curve_fit
from sklearn.preprocessing import minmax_scale
import plotly.graph_objs as go
import pandas as pd
import numpy as np
import ssbtoolkit.Utils as utils

pathways_path=(os.path.join(os.path.split(os.path.realpath(__file__))[0], 'pathways'))

######MAIN CODE#####
class Binding:
    """This class simulate ligand-target binding curves."""
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

class Simulation:
    """
    This class simulates the mathematical models of the signaling pathways.
    """
    class Activation:
        """
        Simulation of the activation of signaling pathways (i.e. activation by agonists)
        """
        def __init__(self):
            self._ligands=None
            self._affinities=None
            self._pathway=None 
            self._receptor_conc=None 
            self._lig_conc_range=None 
            self._ttotal=None 
            self._nsteps=None
            self._binding_kinetics=True 
            self._binding_kinetic_parameters=None
            self.simulation_data=None
            self.processed_data=None

        def SetSimulationParameters(self, **kwargs):
            """
            :parameter ligands:          Required (kwargs list): list of ligands' names (str)
            :parameter affinities:       Required (kwargs list): list of pKd values (flt)
            :parameter pathway:          Required (kwargs str): name of the pathway ('Gs', 'Gi', 'Gq') 
            :parameter receptor_conc:    Required (kwargs flt): receptors concentration (nM)
            :parameter lig_conc_range:   Required (kwargs array): range of ligands' concentration
            :parameter ttotal:           Required (kwargs int): simulation time (seconds)
            :parameter nsteps:           Required (kwargs int): simulation time step
            :parameter binding_kinetics: Optional (kwargs boolean): default (False)

            
            .. warning:: the order of the lists of ligands names and affinities list must be the same. 
            
            """
            self._ligands= kwargs.pop('ligands')
            if 'affinities' in kwargs:
                self._affinities=kwargs.pop('affinities')
            self._pathway=kwargs.pop('pathway')
            self._receptor_conc=kwargs.pop('receptor_conc') 
            self._lig_conc_range=kwargs.pop('lig_conc_range') 
            self._ttotal=kwargs.pop('ttotal') 
            self._nsteps=kwargs.pop('nsteps')
            self._binding_kinetics=kwargs.pop('binding_kinetics') 
            if 'binding_kinetic_parameters' in kwargs:self._binding_kinetic_parameters=kwargs.pop('binding_kinetic_parameters')
            self._DefaultPathwayParametersDataFrame=pd.DataFrame()

            return 

        def PathwayParameters(self):
            """
            Display table with default pathway parameters.

            .. warning:: this functions requires the qgrid library. It doens't work on Google Colab.
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))

            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def UserPathwayParameters(self, path):
            """
            Import user pathway parameters.

            :parameter path:     Required (kwarg str): directory path
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(path)
            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def PathwayParametersToCSV(self, path):
            """
            Export pathway parameters into CSV format.

            :parameter path:     Required (kwarg str): directory path
            """
            self._DefaultPathwayParametersTable.get_changed_df().to_csv(path, index=False)
            print('saved in:', path)
            return 

        def Reactions(self):
            """
            Display pathway reactions.
            """
            from IPython.display import display, HTML
            display(HTML("<style>.container {width:90% !important}</style>"))
            return pd.read_csv(pathways_path+'/{}_reactions.csv'.format(self._pathway))

        def Run(self):
            '''
            This function runs the pathway simulation and returns the raw simulation data.
            '''

            #Check inputs
            if self._ligands==None: raise TypeError("ligands list undefined.")
            elif self._pathway==None: raise TypeError("pathway name undefined.")
            elif self._binding_kinetics==False and self._affinities==None: raise TypeError("affinity_values_dict undefined.")
            elif self._binding_kinetics==True and self._affinities==None: pass
            elif self._lig_conc_range.any() == False: raise TypeError("lig_conc_range undefined.")
            elif self._ttotal==None: raise TypeError("ttotal undefined.")
            elif self._nsteps==None: raise TypeError("nsteps undefined.")
            elif self._receptor_conc==None: raise TypeError("receptor_conc undefined.")
            else: pass
            
            
            #Check Pathway availability and import it
            available_pathways = ['Gs', 'Gi', 'Gq']
            if self._pathway == 'Gz(Gi)': self._pathway = 'Gi'
            if self._pathway not in available_pathways: raise Exception('Unvailable Pathway. Please, introduce it manually. Pathways available: "Gs", "Gi", "Gq".')
            mypathway = importlib.import_module('.'+self._pathway, package='ssbtoolkit.pathways')
            
            #Get default pathway parameters            
            if  self._DefaultPathwayParametersDataFrame.empty and self._binding_kinetic_parameters==None:
                self._DefaultPathwayParametersDataFrame = pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))
                self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()
            
            elif self._DefaultPathwayParametersDataFrame.empty is False and self._binding_kinetic_parameters is None:
                try: 
                    #extract data from qgrid
                    newparameters = self._DefaultPathwayParametersTable.get_changed_df()
                    self._PathwayParameters = newparameters.set_index('Parameter').iloc[:,0].to_dict()
                except:
                    self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()
            
            elif self._DefaultPathwayParametersDataFrame.empty and self._binding_kinetic_parameters is not None:
                self._DefaultPathwayParametersDataFrame = pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))

            elif self._DefaultPathwayParametersDataFrame.empty is False and self._binding_kinetic_parameters is not None:
                try: 
                    #extract data from qgrid
                    newparameters = self._DefaultPathwayParametersTable.get_changed_df()
                    self._PathwayParameters = {**newparameters.set_index('Parameter').iloc[:,0].to_dict(), **self._binding_kinetic_parameters}
                except:
                    self._PathwayParameters = {**self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict(), **self._binding_kinetic_parameters}

            #Input
            t = pl.geomspace(0.00001, self._ttotal, num=self._nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points
            

            #Output
            simulation_data={}

            #Function
            for ligand in self._ligands:
                ligand_name = os.path.splitext(str(ligand))[0]
                data=[]
                utils.PrintProgressBar(0, len(self._lig_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                ###DANGER ZONE###
                if  self._binding_kinetic_parameters is not None: 
                    self._PathwayParameters = {**self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict(), **self._binding_kinetic_parameters[self._ligands.index(ligand)]}
                ######################

                for idx in range(len(self._lig_conc_range)):

                    ligand_conc = self._lig_conc_range[idx]
                    if self._binding_kinetics == False:
                        #get LR conc
                        parameters = {**self._PathwayParameters, 'R_init':self._receptor_conc}
                        LR_conc_init = utils.CalcOccupancy(self._receptor_conc, ligand_conc, 0, self._affinities[self._ligands.index(ligand)], 0)
                        mymodel = mypathway.network(LR=LR_conc_init, kinetics=False, **parameters)
                        simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                        yout = simres.all
                    
                    elif self._binding_kinetics == True:
                        parameters={**self._PathwayParameters,'R_init':self._receptor_conc, 'L_init':self._lig_conc_range[idx] }
                        mymodel = mypathway.network(kinetics=True, **parameters)
                        simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                        yout = simres.all
                        
                    d1={'ligand_conc':ligand_conc, 'time':t }

                    for idx2 in range(len(mypathway.list_of_observables)):
                        d2={mypathway.list_of_observables[idx2]:yout[mypathway.list_of_observables[idx2]]}
                        d1.update(d2)
                    data.append(d1)
                    utils.PrintProgressBar(idx + 1, len(self._lig_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)
                

                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':ligand_name,}
            self.simulation_data = simulation_data
            return

        def Analysis(self):
            '''
            This function calculates the dose-response effect.
            
            :return: instance of processed_data
            '''
            
            if self.simulation_data == None: raise TypeError('There is no simulation data. simulation.activation.run() must be run first.')
          
            # Define all the lists and dictionaries used in this function
            raw_data=[]
            normalized_data=[]
            dose={}

            #defining concentration range
            lig_conc_min = self._lig_conc_range.min()
            lig_conc_max = self._lig_conc_range.max()
            
            #Main function
            for ligand in self.simulation_data:

                #definig and dictionaries used in this loop:
                raw_data_dict={}
                normalized_data_dict={}
            
                # Calculate dose-response curve
                #get metabolite concentration, rescale, and transform data if pathway/metabolite decrease
                # metabolite_raw is not normalized
                if self._pathway == 'Gi' or self._pathway == 'Gz(Gi)':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(1-np.array(metabolite_conc_raw))
                elif self._pathway == 'Gs':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                elif self._pathway == 'Gq':
                    metabolite='IP3'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                else: raise Exception('Unvailable Pathway. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')


                ## save results
                raw_data_dict['x']=self._lig_conc_range
                raw_data_dict['y']=metabolite_conc_raw
                raw_data_dict['label']=self.simulation_data[ligand]['label']

                normalized_data_dict['x']=self._lig_conc_range
                normalized_data_dict['y']=metabolite_conc_norm
                normalized_data_dict['label']=self.simulation_data[ligand]['label']

                ## create a list of all data
                raw_data.append(raw_data_dict)
                normalized_data.append(normalized_data_dict)

                ##Fitting curve to the data

                def equation_dose(X, Bottom, Top, EC50, p):
                    return Bottom + (Top-Bottom)/(1+np.power((EC50/X),p))

                popt_EC50, pcov = curve_fit(equation_dose, self._lig_conc_range, metabolite_conc_norm, bounds=([np.min(metabolite_conc_norm),-np.inf,-np.inf, 0.5],[np.inf,np.max(metabolite_conc_norm),np.inf, 2.5]))

                xfit_EC50 = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
                yfit_EC50 = equation_dose(xfit_EC50, *popt_EC50)

                fit_EC50={'x':xfit_EC50, 'y':yfit_EC50, 'label':self.simulation_data[ligand]['label']}

                dose[ligand] = {'raw_data': raw_data_dict,
                                'normalized_data':normalized_data_dict ,
                                'fitted_data': fit_EC50,
                                'EC50 (μM)': round(popt_EC50[2],5),
                                'pEC50': round(-np.log10(popt_EC50[2]*1E-6),2)}
                
            self.processed_data=dose
            return 

        def ShowCurve(self, save=False, filename=None):
            '''
            Plots the dose-response curve.
            '''

            if self.simulation_data == None: raise TypeError('There is no simulation data. simulation.activation.run() must be run first.')

            import plotly
            import plotly.graph_objs as go
            import plotly.offline as pyoff

            colors = plotly.colors.DEFAULT_PLOTLY_COLORS

            plot_data=[]

            color_id=0
            for ligand in self.processed_data:
                trace_norm = go.Scatter(x=self.processed_data[ligand]['normalized_data']['x'],
                                        y=minmax_scale(self.processed_data[ligand]['normalized_data']['y'])*100 ,
                                        mode='markers',
                                        showlegend=True,
                                        name=self.processed_data[ligand]['normalized_data']['label'],
                                        marker=dict(color=colors[color_id]))
                plot_data.append(trace_norm)

                trace_fitted = go.Scatter(x=self.processed_data[ligand]['fitted_data']['x'],
                                    y=minmax_scale(self.processed_data[ligand]['fitted_data']['y'])*100,
                                    mode='lines',
                                    showlegend=False,
                                    name=self.processed_data[ligand]['fitted_data']['label'],
                                    line=dict(color=colors[color_id]))
                plot_data.append(trace_fitted)
                color_id +=1

            layout = dict(title = '',
                            xaxis = dict(
                                title = '[ligand] μM',
                                type ='log',
                                #range = [-3, 2],
                                exponentformat='e',
                                titlefont=dict(
                                    size=20
                                ),
                                tickfont=dict(
                                    size=20
                                )),
                            yaxis = dict(
                                title = '% Response',
                                #range = [0, 100],
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

            fig = go.Figure(data=plot_data, layout=layout)
            #fig['layout']['yaxis'].update(autorange = True)

            if save==True:
                if filename==None: 
                    filename='plot.html'
                    return pyoff.plot(fig, filename=filename)
                else:
                    ext = os.path.splitext(filename)[-1]
                    if ext == '.png': fig.write_image(filename, scale=3)
                    elif ext == '.html': pyoff.plot(fig, filename=filename)
                    else: raise TypeError("extension not valid. Use png or html.")
            elif save ==False: return fig
            return 
            
        def ShowPotency(self):
            '''
            Return the potency values as a pandas DataFrame.
            '''
            import pandas as pd
            data = Simulation.Activation.PotencyToDict(self)
            df = pd.DataFrame.from_dict(data, orient='index')
            return df

        def PotencyToDict(self):
            '''
            Convert potencies into a dictionary.
            '''
            
            #dependencies
            if self.processed_data == None: raise TypeError('Simulation data unprocessed. simulation.activation.analysis() must be run first.')

            kvalues={}
            for ligand in self.processed_data:
                IC50 = list(self.processed_data[ligand].keys())[-2]
                IC50_value = self.processed_data[ligand][IC50]
                pIC50 = list(self.processed_data[ligand].keys())[-1]
                pIC50_value = self.processed_data[ligand][pIC50]
                kvalues[ligand]={IC50:IC50_value, pIC50:pIC50_value}
            return kvalues
           
        def PotencyToCSV(self, path):
            '''
            Exports the potency values into csv format.

            :parameter path: Required (kwarg str): directory path to save the csv file
            '''

            data = Simulation.Activation.PotencyToDict(self)
            df = pd.DataFrame.from_dict(data, orient='index')
            df.to_csv(path, index=False)
            return
   
    class Inhibition:
        """
        Simulation of the inhibition of signaling pathways (i.e. inhibition by antagonists).
        """
        def __init__(self):
            self._agonist=None
            self._agonist_affinity=None
            self._agonist_submaximal_conc=None
            self._antagonists=None
            self._antagonists_affinities=None
            self._pathway=None 
            self._receptor_conc=None 
            self._lig_conc_range=None 
            self._ttotal=None 
            self._nsteps=None
            self._binding_kinetics=False 
            self._binding_kinetic_parameters=None
            self.simulation_data=None
            self.processed_data=None

        def SetSimulationParameters(self, **kwargs):
            """
            :parameter agonist:                 Required (kwargs str): agonist name 
            :parameter agonist_affinity:        Required (kwargs flt): agonist pKd value
            :parameter agonist_submaximal_conc: Required (kwargs flt): agonist submaximal concentration
            :parameter antagonists:             Required (kwargs list):list of antagonists names (str) 
            :parameter antagonists_affinities:  Required (kwargs list): list of antagonists affinity values (flt)
            :parameter antagonists_conc_range:  Required (kwargs array): range of ligands' concentration (nM)
            :parameter pathway:                 Required (kwargs str): name of the pathway ('Gs', 'Gi', 'Gq')
            :parameter receptor_conc:           Required (kwargs flt): receptors concentration (nM)
            :parameter ttotal:                  Required (kwargs int): simulation time (seconds)
            :parameter nsteps:                  Required (kwargs int): simulation time step
            :parameter kinetics:                Optional (kwargs boolean): default (False)


            :return: instances of all parameters
            
            .. warning:: the order of the lists of the antagonists names and affinities list must be the same. 
            
            """
            self._agonist= kwargs.pop('agonist')
            self._agonist_affinity=kwargs.pop('agonist_affinity')
            self._agonist_submaximal_conc=kwargs.pop('agonist_submaximal_conc')
            self._antagonists=kwargs.pop('antagonists')
            self._antagonists_affinities=kwargs.pop('antagonists_affinities')
            self._pathway=kwargs.pop('pathway')
            self._receptor_conc=kwargs.pop('receptor_conc') 
            self._lig_conc_range=kwargs.pop('lig_conc_range') 
            self._ttotal=kwargs.pop('ttotal') 
            self._nsteps=kwargs.pop('nsteps')
            if 'kinetics' in kwargs: 
                self._binding_kinetics=kwargs.pop('kinetics')
                if self._binding_kinetics==True: raise TypeError("The of Kinetic parameters during an inhibition simulation it is not supported yet.") 
            else: self._binding_kinetics=False
            if 'binding_kinetic_parameters' in kwargs:
                self._binding_kinetic_parameters=kwargs.pop('binding_kinetic_parameters')
            self._DefaultPathwayParametersDataFrame=pd.DataFrame()

            return 

        def PathwayParameters(self):
            """
            Display table with default pathway parameters.

            .. warning:: this functions requires the qgrid library. It doens't work on Google Colab.
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))

            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def UserPathwayParameters(self, path):
            """
            Import user pathway parameters.

            :parameter path:     Required (kwarg str): directory path
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(path)
            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def PathwayParametersToCSV(self, path):
            """
            Export pathway parameters into CSV format.

            :parameter path:     Required (kwarg str): directory path
            """
            self._DefaultPathwayParametersTable.get_changed_df().to_csv(path, index=False)
            print('saved in:', path)
            return 
        
        def Reactions(self):
            """
            Display pathway reactions.
            """
            from IPython.display import display, HTML
            display(HTML("<style>.container {width:90% !important}</style>"))
            return pd.read_csv(pathways_path+'/{}_reactions.csv'.format(self._pathway))

        def Run(self):
            '''
            This function runs the pathway simulation and returns the raw simulation data.
            '''

            #Check inputs
            if self._agonist==None: raise TypeError("agonist undefined.")
            elif self._agonist_affinity==None: raise TypeError("agonist_affinity undifined.")
            elif self._antagonists==None: raise TypeError("antagonists list undefined.")
            elif self._antagonists_affinities==None: raise TypeError("antagonists affinity values undefined.")
            elif self._pathway==None: raise TypeError("pathway undefined.")
            elif self._lig_conc_range.any() == False: raise TypeError("lig_conc_range undefined.")
            elif self._agonist_submaximal_conc == None: raise TypeError("agonist_submaximal_conc undifined.")
            elif self._ttotal==None: raise TypeError("ttotal undefined.")
            elif self._nsteps==None: raise TypeError("nsteps undefined.")
            elif self._receptor_conc==None: raise TypeError("receptor_conc undefined.")
            elif self._binding_kinetics==True: raise TypeError("The of Kinetic parameters during an inhibition simulation it is not supported yet.")
            else: pass

            #check pathway
            available_pathways = ['Gs', 'Gi', 'Gq']
            if self._pathway == 'Gz(Gi)': self._pathway = 'Gi'
            if self._pathway not in available_pathways: raise Exception('Unvailable Pathway. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')
            mypathway = importlib.import_module('.'+self._pathway, package='ssbtoolkit.pathways')
            
            #Get default pathway parameters            
            if  self._DefaultPathwayParametersDataFrame.empty:
                self._DefaultPathwayParametersDataFrame = pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))
                self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()
            
            elif self._DefaultPathwayParametersDataFrame.empty is False:
                try: 
                    #extract data from qgrid
                    newparameters = self._DefaultPathwayParametersTable.get_changed_df()
                    self._PathwayParameters = newparameters.set_index('Parameter').iloc[:,0].to_dict()
                except:
                    self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()            

            #Input
            t = pl.geomspace(0.00001, self._ttotal, num=self._nsteps) # (a,b,c); a is the starting time ; b is the total time simulated ; c is the number of points

            #Output
            simulation_data={}

            #Function
            for ligand in self._antagonists:
                ligand_name = os.path.splitext(ligand)[0]
                data=[]
                utils.PrintProgressBar(0, len(self._lig_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                for idx in range(len(self._lig_conc_range)):

                    ligand_conc = self._lig_conc_range[idx]

                    #get LR conc
                    parameters = {**self._PathwayParameters, 'R_init':self._receptor_conc}
                    LR_conc_init = utils.CalcOccupancy(self._receptor_conc, self._agonist_submaximal_conc, ligand_conc, self._agonist_affinity, self._antagonists_affinities[self._antagonists.index(ligand)])
                    mymodel = mypathway.network(LR=LR_conc_init, kinetics=False, **parameters)
                    simres = ScipyOdeSimulator(mymodel, tspan=t, compiler='cython').run()
                    yout = simres.all

                    d1={'ligand_conc':ligand_conc, 'time':t }

                    for idx2 in range(len(mypathway.list_of_observables)):
                        d2={mypathway.list_of_observables[idx2]:yout[mypathway.list_of_observables[idx2]]}
                        d1.update(d2)
                    data.append(d1)
                    utils.PrintProgressBar(idx + 1, len(self._lig_conc_range), prefix = "{:<15}".format(ligand_name[:15]), suffix = 'Complete', length = 50)

                simulation_data[ligand_name] = {'sim_data':data,
                                            'label':self._agonist+' + ' + ligand_name}

            self.simulation_data=simulation_data
            return

        def Analysis(self):
            '''
            This function calculates the dose-response effect.
            
            :return: instance processed_data
            '''
            #dependencies
            if self.simulation_data == None: raise TypeError('There is no simulation data. simulation.inhibition.run() must be run first.')

            # Define all the lists and dictionaries used in this function
            raw_data=[]
            normalized_data=[]
            dose={}

            #defining concentration range
            #defining concentration range
            lig_conc_min = self._lig_conc_range.min()
            lig_conc_max = self._lig_conc_range.max()

            #Main function
            for ligand in self.simulation_data:

                #definig and dictionaries used in this loop:
                raw_data_dict={}
                normalized_data_dict={}

                # Calculate dose-response curve
                #get metabolite concentration, rescale, and transform data if pathway/metabolite decrease
                # metabolite_raw is not normalized
                if self._pathway == 'Gi' or self._pathway == 'Gz(Gi)':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(1-np.array(metabolite_conc_raw))
                elif self._pathway == 'Gs':
                    metabolite='cAMP'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                elif self._pathway == 'Gq':
                    metabolite='IP3'
                    metabolite_conc_raw=[]
                    for i in range(len(self._lig_conc_range)):
                        n=np.amax(self.simulation_data[ligand]['sim_data'][i]['obs_'+metabolite]) #cad
                        metabolite_conc_raw.append(n)
                    metabolite_conc_norm = minmax_scale(np.array(metabolite_conc_raw))
                else: raise Exception('Unvailable Pathway. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')


                ## save results
                raw_data_dict['x']=self._lig_conc_range
                raw_data_dict['y']=metabolite_conc_raw
                raw_data_dict['label']=self.simulation_data[ligand]['label']

                normalized_data_dict['x']=self._lig_conc_range
                normalized_data_dict['y']=metabolite_conc_norm
                normalized_data_dict['label']=self.simulation_data[ligand]['label']

                ## create a list of all data
                raw_data.append(raw_data_dict)
                normalized_data.append(normalized_data_dict)

                ##Fitting curve to the data

                def equation_dose(X, Bottom, Top, EC50, p):
                    return Bottom + (Top-Bottom)/(1+np.power((EC50/X),p))

                popt_IC50, pcov = curve_fit(equation_dose, self._lig_conc_range, metabolite_conc_norm, bounds=([np.min(metabolite_conc_norm),-np.inf,-np.inf, 0.5],[np.inf,np.max(metabolite_conc_norm),np.inf, 2.5]))

                xfit_IC50 = np.geomspace(lig_conc_min, lig_conc_max, 50000) # These values are the same as the values for the simulation time and not ligand concentration
                yfit_IC50 = equation_dose(xfit_IC50, *popt_IC50)

                fit_IC50={'x':xfit_IC50, 'y':yfit_IC50, 'label':self.simulation_data[ligand]['label']}

                dose[ligand] = {'raw_data': raw_data_dict,
                                'normalized_data':normalized_data_dict ,
                                'fitted_data': fit_IC50,
                                'IC50 (μM)': round(popt_IC50[2],5),
                                'pIC50': round(-np.log10(popt_IC50[2]*1E-6),2)}

            self.processed_data=dose
            return 

        def ShowCurve(self, save=False, filename=None):
            '''
            Plot the dose-response curve.
            '''
            #dependencies
            if self.processed_data == None: raise TypeError('Simulation data unprocessed. simulation.inhibition.analysis() must be run first.')

            import plotly
            import plotly.graph_objs as go
            import plotly.offline as pyoff

            colors = plotly.colors.DEFAULT_PLOTLY_COLORS

            plot_data=[]

            color_id=0
            for ligand in self.processed_data:
                trace_norm = go.Scatter(x=self.processed_data[ligand]['normalized_data']['x'],
                                        y=minmax_scale(self.processed_data[ligand]['normalized_data']['y'])*100 ,
                                        mode='markers',
                                        showlegend=True,
                                        name=self.processed_data[ligand]['normalized_data']['label'],
                                        marker=dict(color=colors[color_id]))
                plot_data.append(trace_norm)

                trace_fitted = go.Scatter(x=self.processed_data[ligand]['fitted_data']['x'],
                                    y=minmax_scale(self.processed_data[ligand]['fitted_data']['y'])*100,
                                    mode='lines',
                                    showlegend=False,
                                    name=self.processed_data[ligand]['fitted_data']['label'],
                                    line=dict(color=colors[color_id]))
                plot_data.append(trace_fitted)
                color_id +=1

            layout = dict(title = '',
                            xaxis = dict(
                                title = '[ligand] μM',
                                type ='log',
                                #range = [-4, 2],
                                exponentformat='e',
                                titlefont=dict(
                                    size=20
                                ),
                                tickfont=dict(
                                    size=20
                                )),
                            yaxis = dict(
                                title = '% Response',
                                #range = [0, 100],
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

            fig = go.Figure(data=plot_data, layout=layout)
            if save==True:
                if filename==None: 
                    filename='plot.html'
                    return pyoff.plot(fig, filename=filename)
                else:
                    ext = os.path.splitext(filename)[-1]
                    if ext == '.png': fig.write_image(filename, scale=3)
                    elif ext == '.html': pyoff.plot(fig, filename=filename)
                    else: raise TypeError("extension not valid. Use png or html.")
            elif save ==False: return fig
            return 

        def ShowPotency(self):
            '''
            Return the potency values as a pandas DataFrame.
            '''
            import pandas as pd
            data = Simulation.Inhibition.PotencyToDict(self)
            df = pd.DataFrame.from_dict(data, orient='index')
            return df

        def PotencyToDict(self):
            '''
            Convert potencies into a dictionary.
            '''
            #dependencies
            if self.processed_data == None: raise TypeError('Simulation data unprocessed. simulation.inhibition.analysis() must be run first.')

            kvalues={}
            for ligand in self.processed_data:
                IC50 = list(self.processed_data[ligand].keys())[-2]
                IC50_value = self.processed_data[ligand][IC50]
                pIC50 = list(self.processed_data[ligand].keys())[-1]
                pIC50_value = self.processed_data[ligand][pIC50]
                kvalues[ligand]={IC50:IC50_value, pIC50:pIC50_value}
            return kvalues
        
        def PotencyToCSV(self, path):
            '''
            Exports the potency values into csv format.

            :parameter path: Required (kwarg str): directory path to save the csv file
            '''
            data = Simulation.Inhibition.PotencyToDict(self)
            df = pd.DataFrame.from_dict(data, orient='index')
            df.to_csv(path, index=False)
            return

    class FitModel:
        """
        Fit a model to experimental data.

        .. note:: This class was developed to reproduce data from a specific experimental setup. Please see tutorial 4 (OXTR pathay). Use carefully!
        """
        def __init__(self):
        
            #fitting parameters
            self._expratio = None
            self._seed = None
            self._maxiter = None
            self._seed_incrementor = None
            self._target_parameter = None
            
            #Pathway parameters
            self._ttotal = None
            self._nsteps = None
            self._pathway = None
            self._observable = None
            self.pathway_parameters = {}
            
        def SetSimulationParameters(self, **kwargs):
            """
            :parameter pathway_parameters: Required (kwargs): dict of pathway parameters
            :parameter pathway:            Required (kwargs str): name of the pathway ('Gs', 'Gi', 'Gq')
            :parameter ttotal:             Required (kwargs int): simulation time (seconds)
            :parameter nsteps:             Required (kwargs int): simulation time step
            :parameter observable:         Required (kwargs str): molecular specie to be measured

            :return: instances of all parameters
                        
            """

            if 'pathway_parameters' in kwargs: 
                self.pathway_parameters = kwargs.pop('pathway_parameters')
                #print('pathway_parameters YES')
            self._DefaultPathwayParametersDataFrame=pd.DataFrame()
            if 'ttotal' in kwargs: 
                self._ttotal = int(kwargs.pop('ttotal'))
                print('ttotal =', self._ttotal)
            else: raise TypeError("ttotal undefined.")
                
            if 'nsteps' in kwargs: 
                self._nsteps = int(kwargs.pop('nsteps', 1000))
                print('nsteps =', self._nsteps)
            
            if 'pathway' in kwargs: 
                self._pathway = str(kwargs.pop('pathway'))
                print('pathway ->', self._pathway)
            else: raise TypeError("pathway undefined.")

            available_pathways = ['Gs', 'Gi', 'Gq', 'OXTR_pathway']
            if self._pathway == 'Gz(Gi)': self._pathway = 'Gi'
            if self._pathway not in available_pathways: raise Exception('Unvailable Pathway. Please, introduce it manually. Networs available: "Gs", "Gi", "Gq".')

            
            if'observable' in kwargs: 
                self._observable = str(kwargs.pop('observable'))
                print('observable ->', self._observable)
            else: raise TypeError("observable undefined.")
            
            return

        def PathwayParameters(self):
            """
            Display table with default pathway parameters.

            .. warning:: this functions requires the qgrid library. It doens't work on Google Colab.
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))

            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def UserPathwayParameters(self, path):
            """
            Import user pathway parameters.

            :parameter path:     Required (kwarg str): directory path
            """
            import qgrid
            self._DefaultPathwayParametersDataFrame =  pd.read_csv(path)
            col_opts = { 'editable': False, 'sortable':False}
            col_defs = {'Value': { 'editable': True, 'width': 150 }}
            self._DefaultPathwayParametersTable = qgrid.show_grid(self._DefaultPathwayParametersDataFrame, column_options=col_opts,column_definitions=col_defs)
            return self._DefaultPathwayParametersTable

        def PathwayParametersToCSV(self, path):
            """
            Export pathway parameters into CSV format.

            :parameter path:     Required (kwarg str): directory path
            """
            self._DefaultPathwayParametersTable.get_changed_df().to_csv(path, index=False)
            print('saved in:', path)
            return 

        def Reactions(self):
            """
            Display pathway reactions.
            """
            from IPython.display import display, HTML
            display(HTML("<style>.container {width:90% !important}</style>"))
            return pd.read_csv(pathways_path+'/{}_reactions.csv'.format(self._pathway))
          
        def Run(self, **kwargs):
            """
            Fits of the model to experimental data.
            
            :parameter expratio:         Required (kwargs flt): experimental signalling specie concentration ratio
            :parameter target_parameter: Required (kwargs str):kinetic parameter to me modified
            :parameter maxiter:          Required (kwargs int): maximum number of iteration
            :parameter seed:             Required (kwargs flt): ramdom seed for scaling the modified parameter
            :parameter seed_incrementor: Required (kwargs flt): seed incrementor (each iteration will increment the seed by this value)
            :parameter seed_decrementor: Required (kwargs flt): seed decrementor (each iteration will decrement the seed by this value)
                        
            """

            from scipy.signal import find_peaks
            import decimal
            #fitting parameters
            if 'expratio' in kwargs: 
                self._expratio = float(kwargs.pop('expratio'))
                print('expratio =', self._expratio)
            else: raise TypeError("exratio undefined.")
                
            if 'seed' in kwargs: 
                self._seed = float(kwargs.pop('seed'))
                print('seed =', self._seed)
            else: raise TypeError("seed undefined.")
                
            if 'maxiter' in kwargs: 
                self._maxiter = int(kwargs.pop('maxiter', 100))
                print('maxiter =', self._maxiter)
    
            if 'seed_incrementor' in kwargs: 
                self._seed_incrementor = float(kwargs.pop('seed_incrementor', 0.1))
                print('seed_incrementor =', self._seed_incrementor)

            if 'seed_decrementor' in kwargs: 
                self._seed_decrementor = float(kwargs.pop('seed_decrementor', 0.1))
                print('seed_decrementor =', self._seed_decrementor)
                
            if 'target_parameter' in kwargs:
                self._target_parameter = str(kwargs.pop('target_parameter'))
                print('target_parameter ->', self._target_parameter)
            else: raise TypeError("target_parameter undefined.")
            
            
            #Get default pathway parameters            
            if  self._DefaultPathwayParametersDataFrame.empty and self.pathway_parameters==None:
                self._DefaultPathwayParametersDataFrame = pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))
                self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()    

            elif self._DefaultPathwayParametersDataFrame.empty is False and self.pathway_parameters is None:
                try: 
                    #extract data from qgrid
                    newparameters = self._DefaultPathwayParametersTable.get_changed_df()
                    self._PathwayParameters = newparameters.set_index('Parameter').iloc[:,0].to_dict()
                except:
                    self._PathwayParameters = self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict()
            
            elif self._DefaultPathwayParametersDataFrame.empty and self.pathway_parameters is not None: 
                self._DefaultPathwayParametersDataFrame = pd.read_csv(pathways_path+'/{}_parameters.csv'.format(self._pathway))
                self._PathwayParameters = {**self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict(), **self.pathway_parameters}

            elif self._DefaultPathwayParametersDataFrame.empty is False and self.pathway_parameters is not None:
                try: 
                    #extract data from qgrid
                    newparameters = self._DefaultPathwayParametersTable.get_changed_df()
                    self._PathwayParameters = {**newparameters.set_index('Parameter').iloc[:,0].to_dict(), **self.pathway_parameters}
                except:
                    self._PathwayParameters = {**self._DefaultPathwayParametersDataFrame.set_index('Parameter').iloc[:,0].to_dict(), **self.pathway_parameters}



            #simulation parameters:
            if not self._ttotal: 
                raise TypeError("simulation parameters unknown. Set the the simulation parameters first wiht set_simulation_parameters()")
            
            #Main function
            mypathway = importlib.import_module('.'+self._pathway, package='ssbtoolkit.pathways')
            self.simtime = pl.geomspace(0.00001, self._ttotal, num=self._nsteps) 

            #Simulation 1
            pathway_model = mypathway.network(LR=None, kinetics=True, **self._PathwayParameters)
            sim1 = ScipyOdeSimulator(pathway_model, tspan=self.simtime,compiler='cython').run()
            self.simres1 = sim1.all
                                
            def calc_ratio(self):
                
                
                #Simulation 2
                sim2 = ScipyOdeSimulator(mypathway.network(**self.new_pathway_parameters), tspan=self.simtime, compiler='cython').run()
                self.simres2 = sim2.all

                #analysis
                obs_name = 'obs_'+self._observable
                obs_1  = self.simres1[obs_name]
                obs_2  = self.simres2[obs_name]
                
                if 'time_in' in self.new_pathway_parameters:
                    self._time = np.take(self.simtime, np.where(self.simtime > int(self.new_pathway_parameters['time_in'])))[0]
                    obs_curve_1  = np.take(obs_1,  np.where(self.simtime > int(self.new_pathway_parameters['time_in'])))[0]
                    obs_curve_2  = np.take(obs_2, np.where(self.simtime > int(self.new_pathway_parameters['time_in'])))[0]
                                
                else:
                    self._time = np.take(self.simtime, np.where(self.simtime>0))[0]
                    obs_curve_1  = np.take(obs_1,  np.where(self.simtime > 0))[0]
                    obs_curve_2  = np.take(obs_2, np.where(self.simtime > 0))[0]
                        
                obs_peaks_1, _  = find_peaks(obs_curve_1)
                obs_peaks_2, _  = find_peaks(obs_curve_2)

                vmax_obs_curve_1  = obs_curve_1[obs_peaks_1][-1]*1E3
                vmax_obs_curve_2  = obs_curve_2[obs_peaks_2][-1]*1E3

                obs_ratio = round(vmax_obs_curve_2/vmax_obs_curve_1, abs(decimal.Decimal(str(self._expratio).rstrip('0')).as_tuple().exponent))

                self._obs_curve_1=obs_curve_1
                self._obs_curve_2=obs_curve_2
                self._obs_peaks_1=obs_peaks_1
                self._obs_peaks_2=obs_peaks_2
                self._vmax_obs_curve_1=vmax_obs_curve_1
                self._vmax_obs_curve_2=vmax_obs_curve_2

                return obs_ratio
        
            self._iteration=1
            print('\n')
            
            self._lst_ratio=[]
            self._lst_seed=[]
        
            for idx in range(self._maxiter):
                
                prefix = 'iteration'
                iteration_n = str(self._iteration)
                print(f'\r{prefix} {iteration_n}', end='\r')
                
                self.new_pathway_parameters={**self._PathwayParameters, **{self._target_parameter:mypathway.defaultParameters[self._target_parameter]*self._seed}}
                self.obs_ratio = calc_ratio(self)

                if self.obs_ratio == self._expratio:
                    self._lst_ratio.append(self.obs_ratio)
                    self._lst_seed.append(self._seed)
                    self._fold=round(self._seed, abs(decimal.Decimal(str(self._expratio).rstrip('0')).as_tuple().exponent))
                    print('\n\nDONE!\n', '\nRatio: '+str(self.obs_ratio), '\nFOLD: '+str(self._fold), '\nNumber of iterations: '+str(self._iteration))
                    break
                elif self.obs_ratio < self._expratio:
                    
                    self._lst_ratio.append(self.obs_ratio)
                    self._lst_seed.append(self._seed)
                    self._iteration+=1
                    self._seed += self._seed_incrementor     

                else:
                    self._lst_ratio.append(self.obs_ratio)
                    self._lst_seed.append(self._seed)


                    self._iteration+=1
                    self._seed -= self._seed_decrementor
    
            return

        def PlotIterations(self, save=False, filename=None):
            '''
            Plot iterations. 
            '''
            import plotly.offline as pyoff
            #dependencies
            if self._iteration == None: raise TypeError('Simulation data not exist. simulation.fitModel.run() must be run first.')

            #import plotly
            import plotly.graph_objs as go

            iterations = np.arange(1,self._iteration+1)

            trace=dict(type='scatter', x=self._lst_seed, y=self._lst_ratio, mode='markers', 
                    marker=dict(color= iterations, colorscale='Bluered_r', size=14, colorbar=dict(thickness=20, title='iteration number')))
            #axis_style=dict(zeroline=False, showline=True, mirror=True)
            layout = dict(title = '',
                xaxis = dict(
                    title = 'seed',
                    titlefont=dict(
                        size=20
                    ),
                    tickfont=dict(
                        size=20
                    )),
                yaxis = dict(
                    title = '['+self._observable+']' + ' ratio',
                    titlefont=dict(
                        size=20),
                    tickfont=dict(
                    size=20)

                ),
                legend=dict(font=dict(size=15)),
                autosize=False,
                width=850,
                height=650
                )

            fig = go.Figure(data=[trace], layout=layout)
            if save==True:
                if filename==None: 
                    filename='plot.html'
                    return pyoff.plot(fig, filename=filename)
                else:
                    ext = os.path.splitext(filename)[-1]
                    if ext == '.png': fig.write_image(filename, scale=3)
                    elif ext == '.html': pyoff.plot(fig, filename=filename)
                    else: raise TypeError("extension not valid. Use png or html.")
            elif save ==False: return fig
            return fig

        def ShowGraphs(self, save=False, filename=None):
            '''
            Plot the amount of obeservable in function of time, Amplitude, Area Under the Curve, and Full Width at Half Maximum. 

            :parameter save:     Optional (kwarg boolean): default False
            :parameter filename: Optional (kwarg str)
            '''

            from IPython.core.display import display, HTML
            display(HTML("<style>.container { width:90% !important; }</style>"))


            from plotly.subplots import make_subplots
            from scipy.signal import peak_widths
            from sklearn import metrics
            import plotly.offline as pyoff


            half_1 = peak_widths(self._obs_curve_1, self._obs_peaks_1, rel_height=0.5)
            half_2 = peak_widths(self._obs_curve_2, self._obs_peaks_2, rel_height=0.5)
            fwhm_1 = self._time[int(half_1[3])]-self._time[int(half_1[2])]
            fwhm_2 = self._time[int(half_2[3])]-self._time[int(half_2[2])]


            fig = make_subplots(rows=2, cols=2,vertical_spacing=0.15,
                                subplot_titles=("{} concentration".format(self._observable), "Amplitude", "Area under the curve", "Full Width at Half Maximum"))

            ####################
            #### MAIN PLOT  ####
            ####################
            fig.add_trace(go.Scatter(x=self._time, y=self._obs_curve_1*1E3, name='control'), row=1, col=1)
            fig.add_trace(go.Scatter(x=self._time, y=self._obs_curve_2*1E3, name='{}-fold'.format(self._fold)), row=1, col=1)
            fig.add_trace(go.Scatter(x=self._time[self._obs_peaks_1], y=self._obs_curve_1[self._obs_peaks_1]*1E3,
                                    name='max value', showlegend=False, mode='markers', 
                                    marker=dict(symbol='x', size=13, color='Black')), row=1,col=1)
            fig.add_trace(go.Scatter(x=self._time[self._obs_peaks_2], y=self._obs_curve_2[self._obs_peaks_2]*1E3,
                                    name='max value', showlegend=False, mode='markers', 
                                    marker=dict(symbol='x', size=13, color='Black')), row=1,col=1)
            fig.add_shape(type='line', x0=self._time[int(half_1[2])],y0=half_1[1][0]*1E3, x1=self._time[int(half_1[3])], y1=half_1[1][0]*1E3,
                        line=dict(color='Blue',dash='dash'),xref='x',yref='y', row=1, col=1)
            fig.add_shape(type='line', x0=self._time[int(half_2[2])],y0=half_2[1][0]*1E3, x1=self._time[int(half_2[3])], y1=half_2[1][0]*1E3,
                        line=dict(color='Red',dash='dash'),xref='x',yref='y', row=1, col=1)

            # Update xaxis properties
            fig.update_xaxes(title_text="Time (s)", showgrid=False, row=1, col=1, titlefont=dict(size=18), 
                            linecolor='black', linewidth=2,
                            ticks='inside', tickfont=dict(size=18), tickcolor='black', ticklen=10, tickwidth=2)

            fig.update_yaxes(title_text=self._observable+' (nM)', titlefont=dict(size=18), showgrid=False, row=1, col=1, 
                            linecolor='black', linewidth=2,
                            ticks='inside', tickfont=dict(size=18), tickcolor='black', ticklen=10, tickwidth=2)


            ####################
            #### AMPLITUDE  ####
            ####################

            AMP_labels = [1,2]
            AMP_values = [self._vmax_obs_curve_1, self._vmax_obs_curve_2]
            fig.add_trace(go.Bar(x=AMP_labels,y=AMP_values, width = [0.35,0.35], showlegend=False, marker_color='black', name=''), row=1, col=2 )

                

            # Update xaxis properties
            fig.update_xaxes(row=1, col=2, showgrid=False, linecolor='black', linewidth=2, range=[0,3],
                            tickmode='array', tickvals=[1,2], ticktext=['control', '{}-fold'.format(self._fold)], tickfont=dict(size=18))

            fig.update_yaxes(showgrid=False, range=[round((min(AMP_values)-min(AMP_values)*0.5)/5)*5,round((max(AMP_values)+max(AMP_values)*0.5)/5)*5 ], row=1, col=2,
                            title_text=self._observable+' (nM)', titlefont=dict(size=18),
                            linecolor='black', linewidth=2, ticks='inside', ticklen=10, tickwidth=2, tickfont=dict(size=18))

            # Add diff lines
            AMP_diffs = [max(AMP_values) - v for v in AMP_values]
            AMP_diff_labels = dict(zip(AMP_labels, AMP_diffs))
            fig.add_trace(go.Scatter(name='',x=[1,1.5,2], y=[max(AMP_values)+(max(AMP_values)*0.3)]*3, mode = 'lines+text',showlegend=False, line=dict(color='black', width=1),text=['', 'diff. = {} nM'.format(round(AMP_diffs[0], 3)),''], textposition='top center'), row=1, col=2)
            fig.add_trace(go.Scatter(name='',x=[AMP_labels[0]-0.175, AMP_labels[0]+0.175], y=[AMP_values[0]+(AMP_values[0]*0.03)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=1, col=2)
            fig.add_trace(go.Scatter(name='',x=[AMP_labels[1]-0.175, AMP_labels[1]+0.175], y=[AMP_values[1]+(AMP_values[1]*0.03)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=1, col=2)
            fig.add_trace(go.Scatter(name='',x=[AMP_labels[0], AMP_labels[0]], y=[AMP_values[0]+(AMP_values[0]*0.03), max(AMP_values)+(max(AMP_values)*0.3)], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=1, col=2)
            fig.add_trace(go.Scatter(name='',x=[AMP_labels[1], AMP_labels[1]], y=[AMP_values[1]+(AMP_values[1]*0.03), max(AMP_values)+(max(AMP_values)*0.3)], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=1, col=2)
            

            ####################
            ####     AUC    ####
            ####################

            # Data
            AUC_labels = [1,2]
            AUC_values = [round(metrics.auc(self._time, self._obs_curve_1),2), round(metrics.auc(self._time, self._obs_curve_2),2)]
            fig.add_trace(go.Bar(x=AUC_labels,y=AUC_values, width = [0.35,0.35], showlegend=False, marker_color='black', name=''), row=2, col=1 )
                    
            # Update xaxis properties
            fig.update_xaxes(row=2, col=1, tickmode='array', showgrid=False, range=[0,3], linecolor='black', linewidth=2,
                            tickvals=[1,2], ticktext=['control', '{}-fold'.format(self._fold)], tickfont=dict(size=18))

            fig.update_yaxes(row=2, col=1,showgrid=False,  title_text=self._observable+' (nM)', range=[round((min(AUC_values)-min(AUC_values)*0.5)/5)*5,round((max(AUC_values)+max(AUC_values)*0.5)/5)*5], 
                            titlefont=dict(size=18),linecolor='black', linewidth=2, 
                            ticks='inside', tickfont=dict(size=18),ticklen=10, tickwidth=2)

            # Add diff lines
            AUC_diffs = [max(AUC_values) - v for v in AUC_values]
            AUC_diff_labels = dict(zip(AUC_labels, AUC_diffs))
            fig.add_trace(go.Scatter(name='',x=[1,1.5,2], y=[max(AUC_values)+(max(AUC_values)*0.3)]*3, mode = 'lines+text',showlegend=False, 
                                    line=dict(color='black', width=1), text=['', 'diff. = {} nM'.format(round(AUC_diffs[0], 3)),''], textposition='top center'), row=2, col=1)
            fig.add_trace(go.Scatter(name='',x=[AUC_labels[0]-0.175, AUC_labels[0]+0.175], y=[AUC_values[0]+(AUC_values[0]*0.03)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=1)
            fig.add_trace(go.Scatter(name='',x=[AUC_labels[1]-0.175, AUC_labels[1]+0.175], y=[AUC_values[1]+(AUC_values[1]*0.03)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=1)
            fig.add_trace(go.Scatter(name='',x=[AUC_labels[0], AUC_labels[0]], y=[AUC_values[0]+(AUC_values[0]*0.03), max(AUC_values)+(max(AUC_values)*0.3)], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=1)
            fig.add_trace(go.Scatter(name='',x=[AUC_labels[1], AUC_labels[1]], y=[AUC_values[1]+(AUC_values[1]*0.03), max(AUC_values)+(max(AUC_values)*0.3)], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=1)


            ####################
            ####    FWHM    ####
            ####################
            # Data
            FWHM_labels = [1,2]
            FWHM_values = [fwhm_1, fwhm_2]
            fig.add_trace(go.Bar(x=FWHM_labels,y=FWHM_values, width = [0.35,0.35], showlegend=False,marker_color='black', name=''), row=2, col=2 )
            
            # Update xaxis properties
            fig.update_xaxes(row=2, col=2, showgrid=False, range=[0,3], linecolor='black', linewidth=2, 
                            tickmode='array', tickvals=[1,2], ticktext=['control', '{}-fold'.format(self._fold)], tickfont=dict(size=18))

            fig.update_yaxes(row=2, col=2, showgrid=False, range=[self.pathway_parameters['time_in'],round((max(FWHM_values)+(max(FWHM_values)-self.pathway_parameters['time_in'])*0.5)/5)*5], 
                            title_text='Time (s)', titlefont=dict(size=18), linecolor='black', linewidth=2,
                            ticks='inside', ticklen=10, tickwidth=2, tickfont=dict(size=18))

            # Add diff lines
            FWHM_diffs = [max(FWHM_values) - v for v in FWHM_values]
            FWHM_diff_labels = dict(zip(FWHM_labels, FWHM_diffs))
            line_height = max(FWHM_values)+((max(FWHM_values)-self.pathway_parameters['time_in'])*0.30)
            fig.add_trace(go.Scatter(x=[1,1.5,2], y=[line_height]*3, mode = 'lines+text',showlegend=False, line=dict(color='black', width=1),name='',
                                    text=['', 'diff. = {} s'.format(round(FWHM_diffs[0], 3)),''], textposition='top center'), row=2, col=2)
            fig.add_trace(go.Scatter(name='',x=[FWHM_labels[0]-0.175, FWHM_labels[0]+0.175], y=[FWHM_values[0]+(FWHM_values[0]*0.005)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=2)
            fig.add_trace(go.Scatter(name='',x=[FWHM_labels[1]-0.175, FWHM_labels[1]+0.175], y=[FWHM_values[1]+(FWHM_values[1]*0.005)]*2, mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=2)
            fig.add_trace(go.Scatter(name='',x=[FWHM_labels[0], FWHM_labels[0]], y=[FWHM_values[0]+(FWHM_values[0]*0.005), line_height], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=2)
            fig.add_trace(go.Scatter(name='',x=[FWHM_labels[1], FWHM_labels[1]], y=[FWHM_values[1]+(FWHM_values[1]*0.005), line_height], mode = 'lines',showlegend=False, line=dict(color='black', width=1)), row=2, col=2)


            ####################
            ####   FIGURE   ####
            ####################

            fig.update_layout(height=1200, width=1300, title_text="", plot_bgcolor='white',showlegend=True, 
                            legend=dict(yanchor="top", x=0.3, y=.99,font=dict(family="sans-serif", size=14,color="black")))
            fig.update_annotations(font_size=20, font_color='black')
            
            if save==True:
                if filename==None: 
                    filename='plot.html'
                    return pyoff.plot(fig, filename=filename)
                else:
                    ext = os.path.splitext(filename)[-1]
                    if ext == '.png': fig.write_image(filename, scale=3)
                    elif ext == '.html': pyoff.plot(fig, filename=filename)
                    else: raise TypeError("extension not valid. Use png or html.")
            elif save ==False: return fig
            
            return
