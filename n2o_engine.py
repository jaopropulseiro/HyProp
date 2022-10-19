import math
import numpy as np
from scipy.integrate import solve_ivp
from CoolProp.CoolProp import PropsSI


class Grain:
    def __init__(self,grain_length,grain_inner_diam,
        grain_outer_diam):

        self.grain_length = grain_length
        self.grain_inner_diam = grain_inner_diam
        self.grain_outer_diam = grain_outer_diam

class Injector:

    def __init__(self,injector_diam,nof_injector_holes,injector_Cd):

        self.injector_diam = injector_diam
        self.nof_injector_holes = nof_injector_holes
        self.injector_area = nof_injector_holes*0.25*math.pi*math.pow(injector_diam,2)
        self.injector_Cd = injector_Cd

class Vent:

    def __init__(self,vent_diam,vent_Cd):
        self.vent_diam = vent_diam
        self.vent_area = 0.25*math.pi*math.pow(vent_diam,2)
        self.vent_Cd = vent_Cd

class Nozzle:

    def __init__(self,nozzle_throat_diam,nozzle_outlet_diam):
        self.nozzle_throat_diam = nozzle_throat_diam
        self.nozzle_outlet_diam = nozzle_outlet_diam
        self.nozzle_throat_area = 0.25*math.pi*math.pow(nozzle_throat_diam,2)

class Engine:

    def __init__(
        self,grain,injector,vent,n2o_mass,chamber_vol,tank_volume,
        nozzle
    ):
        self.grain = grain
        self.injector = injector
        self.vent = vent
        self.n2o_mass = n2o_mass
        self.chamber_vol = chamber_vol
        self.tank_volume = tank_volume
        self.nozzle = nozzle
        

class Combustion:

    def __init__(
        self,gas_molecular_mass,fuel_density,
        isentropic_coeff,combustion_temperature,
        combustion_law_coefs,cstar
    ):
        self.gas_molecular_mass = gas_molecular_mass
        self.gas_constant = 8314/self.gas_molecular_mass
        self.fuel_density = fuel_density
        self.isentropic_coeff=isentropic_coeff
        self.combustion_temperature = combustion_temperature
        self.combustion_law_coefs = combustion_law_coefs
        self.cstar = cstar

class InitialConditions:

    def __init__(self,n2o_temperature,n2o_mass,atmospheric_pressure,grain):
        self.n2o_temperature = n2o_temperature
        self.n2o_mass = n2o_mass
        self.chamber_pressure = atmospheric_pressure
        self.atmospheric_pressure = atmospheric_pressure
        self.grain_inner_radius = grain.grain_inner_diam/2.

class NitrousOxide:

    def __init__(self):
        
        return

    def UpdateProperties(
        tank_volume,initial_conditions=None,
        n2o_quality=None,n2o_s=None,n2o_rho=None
    ):
        n2o_state = {}
        if initial_conditions:
            n2o_state['P_tank'] = PropsSI(
                'P',
                'T',InitialConditions.n2o_temperature,
                'Q',0,'NitrousOxide'
            )
            n2o_state['rhov'] = PropsSI(
                'D',
                'T',InitialConditions.n2o_temperature,
                'Q',1,'NitrousOxide'
            )
            n2o_state['rhol'] = PropsSI(
                'D',
                'T',InitialConditions.n2o_temperature,
                'Q',0,'NitrousOxide'
            )
            n2o_state['sv'] = PropsSI(
                'S',
                'T',InitialConditions.n2o_temperature,
                'Q',1,'NitrousOxide'
            )
            n2o_state['sl'] = PropsSI(
                'S',
                'T',InitialConditions.n2o_temperature,
                'Q',0,'NitrousOxide'
            )
            n2o_state['X'] = (
                n2o_state['rhov']*n2o_state['rhol']*tank_volume-n2o_state['rhov']*InitialConditions.n2o_mass
            )/(InitialConditions.n2o_mass*(n2o_state['rhol']-n2o_state['rhov']))
            total_entropy = InitialConditions.n2o_mass*(
                n2o_state['sl']*(1-n2o_state['X']) + n2o_state['sv']*n2o_state['X']
            )
            return n2o_state,total_entropy
        '''else:
            if n2o_quality<1. and n2o_quality>0.:
                n2o_state['P_tank'] = PropsSI(
                    'P',
                    'D',n2o_rho,
                    'S',n2o_s,'NitrousOxide'
                )
                n2o_state['rhov'] = PropsSI(
                    'D',
                    'T',InitialConditions.n2o_temperature,
                    'Q',1,'NitrousOxide'
                )
                n2o_state['rhol'] = PropsSI(
                    'D',
                    'T',InitialConditions.n2o_temperature,
                    'Q',0,'NitrousOxide'
                )
                n2o_state['sv'] = PropsSI(
                    'S',
                    'T',InitialConditions.n2o_temperature,
                    'Q',1,'NitrousOxide'
                )
                n2o_state['sl'] = PropsSI(
                    'S',
                    'T',InitialConditions.n2o_temperature,
                    'Q',0,'NitrousOxide'
                )
                n2o_state['X'] = (
                    n2o_state['rhov']*n2o_state['rhol']*tank_volume-n2o_state['rhov']*InitialConditions.n2o_mass
                )/(InitialConditions.n2o_mass*(n2o_state['rhol']-n2o_state['rhov']))
                total_entropy = InitialConditions.n2o_mass*(
                    n2o_state['sl']*(1-n2o_state['X']) + n2o_state['sv']*n2o_state['X']
                )
            return n2o_state'''
            
    def GasMassFlow(n2o_state,area,Cd,choked_flow):
        gamma = n2o_state['gamma']
        if choked_flow:
            mdot = Cd*area*math.sqrt(
                gamma*n2o_state['P_tank']*n2o_state['rhov']*\
                math.pow(2/(gamma+1),(gamma+1)/(gamma-1))
            )
        else:
            mdot = Cd*area*math.sqrt(
                2*gamma/(gamma+1)*n2o_state['P_tank']*n2o_state['rhov']*(
                    math.pow(n2o_state['P2']/n2o_state['P_tank'],2/gamma) -\
                    math.pow(n2o_state['P2']/n2o_state['P_tank'],(gamma+1)/gamma)
                )
            )
        return mdot
    
    def LiquidMassFlow(n2o_state,area,Cd):
        kappa = math.sqrt((n2o_state['P1']-n2o_state['P2'])/(n2o_state['Pv1']-n2o_state['P2']))
        mdot_inc = Cd*area*\
            math.sqrt(2*n2o_state['rho1']*(n2o_state['P1']-n2o_state['P2']))
        mdot_HEM = Cd*area*n2o_state['rho2']*\
            math.sqrt(2*(n2o_state['h1']-n2o_state['h2']))
        mdot = ( (1-1/(1+kappa)*mdot_inc + 1/(1+kappa)*mdot_HEM ))
        return mdot

class ODE_System:

    def __init__(self):
        return

    def StateEquations(t,state_vector,n2o_state,engine,combustion,atmospheric_pressure):
        # State vector:
        # 0 - n2o mass through injector
        # 1 - n2o mass through vent
        # 2 - combustion chamber pressure
        # 3 - fuel grain inner radius
        gamma = n2o_state['gamma']
        n2o_quality = n2o_state['X']

        pressure_ratio_threshsold = math.pow(1+0.5*(gamma-1),-gamma/(gamma-1))
        sonic_vent = True if atmospheric_pressure/n2o_state['P_tank'] < pressure_ratio_threshsold else False
        sonic_injector = True if state_vector[2]/n2o_state['P_tank'] < pressure_ratio_threshsold else False

        if n2o_quality<1. and n2o_quality>0.:
            mdot_injector = NitrousOxide.LiquidMassFlow(
                n2o_state,engine.injector.injector_area,
                engine.injector.injector_Cd
            )
            mdot_vent = NitrousOxide.GasMassFlow(
                n2o_state,engine.vent.vent_area,
                engine.vent.vent_Cd,sonic_vent
            )
        else:
            mdot_injector = NitrousOxide.GasMassFlow(
                n2o_state,engine.injector.injector_area,
                engine.injector.injector_Cd,sonic_injector
            )
            mdot_vent = NitrousOxide.GasMassFlow(
                n2o_state,engine.vent.vent_area,
                engine.vent.vent_Cd,sonic_vent
            )

        [a,n]=combustion.combustion_law_coefs
        port_area = math.pi*math.pow(state_vector[3],2)
        n2o_mass_flux = mdot_injector/port_area
        regression_rate = a*math.pow(n2o_mass_flux,n)
        burn_area = 2*math.pi*state_vector[3]*engine.grain.grain_length
        grain_vol = engine.grain.grain_length*0.25*math.pi*(
            math.pow(engine.grain.grain_outer_diam,2)-math.pow(2*state_vector[3],2)
        )
        chamber_volume = engine.chamber_vol - grain_vol
        chamber_pressure_diff = burn_area*regression_rate/chamber_volume*(
            combustion.fuel_density*combustion.gas_constant*combustion.combustion_temperature-state_vector[2]
        ) - state_vector[2]*(
            engine.nozzle.nozzle_throat_area/chamber_volume*math.sqrt(
                gamma*combustion.gas_constant*combustion.combustion_temperature*\
                math.pow(2/(gamma+1),(gamma+1)/(gamma-1))
            )
        ) + combustion.gas_constant*combustion.combustion_temperature/chamber_volume*state_vector[0]
        return np.array([
            mdot_injector,mdot_vent,chamber_pressure_diff,regression_rate
        ])

    def Solution_Loop(
        InitialConditions,engine,combustion,
        timestep,residual_mass=300e-3
    ):
        n2o_state,total_entropy = NitrousOxide.UpdateProperties(True,engine.tank_volume,initial_conditions)
        n2o_mass = InitialConditions.n2o_mass
        injector_evac_mass = 0.
        vent_evac_mass = 0.
        chamber_pressure = InitialConditions.atmospheric_pressure
        grain_inner_radius = InitialConditions.grain_inner_radius
        saturation = True
        t = 0
        step_solution = [injector_evac_mass,vent_evac_mass,chamber_pressure,grain_inner_radius]
        while n2o_mass>residual_mass:
            step_solution = solve_ivp(
                ODE_System.StateEquations,
                [t,t+timestep],
                step_solution,
                args=(n2o_state,engine,combustion,InitialConditions.atmospheric_pressure),
                first_step=timestep,max_step=timestep
            )
            n2o_mass -= (injector_evac_mass + vent_evac_mass)
            n2o_state['rho'] = n2o_mass/engine.tank_volume
            if saturation:
                total_entropy -= (step_solution[0]*n2o_state['sl']+step_solution[1]*n2o_state['sv'])
            else:
                total_entropy -= (step_solution[0]+step_solution[1])*n2o_state['s']
            n2o_state['s'] = total_entropy/n2o_mass
            n2o_state['X'] = PropsSI(
                'Q',
                'D',n2o_state['rho'],
                'S',n2o_state['s'],'NitrousOxide'
            )
            n2o_liquid_mass = n2o_state['X']*n2o_mass
            n2o_gas_mass = (1.-n2o_state['X'])*n2o_mass
            t += timestep
            
        return
