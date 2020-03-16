"""Class for converter models."""




class InverterModels:
    """
    SwitcInverter class.
    
    Attributes:
        ():
        
    """
    
    def ODE_model_single_phase_EMT(self,y,t,signals,grid,sim):
        """ODE model of inverter branch."""
        
        self.ia,dummy = y   # unpack current values of y
        
        Vdc = 100.0 #Get DC link voltage
        control_signal = self.control_signal_calc(signals,t)
        self.vta = self.vt_calc(Vdc,control_signal)
        
        self.va = self.va_calc(t,grid,sim.use_grid)
                               
        dia = (1/self.Lf)*(-self.Rf*self.ia -self.va + self.vta)
                
        result =  [dia,dummy]
        
        return np.array(result)
    
    
    def ODE_model_single_phase_dynamicphasor(self,y,t,signals,grid,sim):
        """Dynamic phasor."""
        
        iaR,iaI = y   # unpack current values of y
        Vdc = 100.0 #Get DC link voltage
        winv = 2*math.pi*60
        
        self.ia = iaR + 1j*iaI
        self.vta = self.half_bridge_phasor(Vdc,1.0+1j*0.0)
        
        diaR = (1/self.Lf)*(-self.Rf*self.ia.real - self.Rload*self.ia.real + self.vta.real) + (winv)*self.ia.imag 
        diaI = (1/self.Lf)*(-self.Rf*self.ia.imag - self.Rload*self.ia.imag + self.vta.imag) - (winv)*self.ia.real  
        result =  [diaR,diaI]
        
        return np.array(result)    
    
    
    def single_phase_half_bridge_switching(self,Vdc,S1):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        S11 = self.calc_primary(S1)
        S12 = self.calc_complimentary(S1)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in ideal half bridge.'

        #print('S11:{},S12:{}'.format(S11,S12))
        Van = (S11 - S12)*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
    
    def single_phase_full_bridge_switching(self,Vdc,S1,S2):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        S11 = self.calc_primary(S1)
        S12 = self.calc_complimentary(S1)
        
        S21 = calc_primary(S2)
        S22 = calc_complimentary(S2)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in full bridge.'
        assert S21+S22 == 1, 'S21 and S22 switches cannot be both ON or OFF at the same time in full bridge.'

        Van = (S11 - S12)*(self.Vdc/2) - (S21 - S22)*(self.Vdc/2)
        
        #print('Van:{}'.format(Van))

        return Van
    
    
    def three_phase_full_bridge_switching(Vdc,S1,S2,S3):
        """Simulates a bridge in inverter"""

        S11 = calc_primary(S1)
        S12 = calc_complimentary(S1)

        S21 = calc_primary(S2)
        S22 = calc_complimentary(S2)

        S31 = calc_primary(S3)
        S32 = calc_complimentary(S3)
        
        assert S11+S12 == 1, 'S11 and S12 switches cannot be both ON or OFF at the same time in ideal half bridge.'
        assert S21+S22 == 1, 'S21 and S22 switches cannot be both ON or OFF at the same time in ideal half bridge.'
        assert S31+S32 == 1, 'S31 and S32 switches cannot be both ON or OFF at the same time in ideal half bridge.'        

        print('S1:{},S2:{},S3:{}'.format(S11,S21,S31))

        Vno =  (self.Vdc/6)*(2*S11+2*S21+2*S31-3)

        Van = (self.Vdc/2)*(S11-S12)-Vno
        Vbn = (self.Vdc/2)*(S21-S22)-Vno
        Vcn = (self.Vdc/2)*(S31-S32)-Vno

        #Van = (2*S11 - S21 - S31)*(Vdc/3)
        #Vbn = (2*S21 - S11 - S31)*(Vdc/3)
        #Vcn = (2*S31 - S21 - S11)*(Vdc/3)
        print('Vno:{},Van+Vbn+Vcn:{}'.format(Vno,Van+Vbn+Vcn))

        print('Van:{},Vbn:{},Vcn:{}'.format(Van,Vbn,Vcn))
        print('Vab:{},Vbc:{},Vca:{}'.format(Van-Vbn,Vbn-Vcn,Vcn-Van))

        return Van,Vbn,Vcn
    
    
    def single_phase_half_bridge_average(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert m>=-1 and m <= 1, 'duty cycle should be between 0 and 1.'

        Van = m*(self.Vdc/2)

        return Van
    
    def single_phase_full_bridge_average(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert m>=-1 and m <= 1, 'duty cycle should be between 0 and 1.'

        Van = m*(self.Vdc)        
        
        return Van
    
    def single_phase_half_bridge_phasor(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert isinstance(m,complex), 'duty cycle should be complex phasor.'

        Van = m*(self.Vdc/2)        
        
        return Van
       
    def single_phase_full_bridge_phasor(self,Vdc,m):
        """Simulates a bridge in inverter"""
        
        self.update_Vdc(Vdc)

        assert isinstance(m,complex), 'duty cycle should be complex phasor.'

        Van = m*(self.Vdc)        
        
        return Van