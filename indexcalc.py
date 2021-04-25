import numpy as np

def mu_per_gy(data):
    """
    Devuelve el índice de complejidad para el plan especificado en data
    """
    MU_beam=[] # UM de cada haz del tratamiento
    d_beam=[]  # Dosis de cada haz del tratamiento

    for fraction in data.FractionGroupSequence:
        for beam in fraction.ReferencedBeamSequence:
            MU_beam.append(beam.BeamMeterset)
            d_beam.append(beam.BeamDose)

    #esto es justamente el ínidce de complejidad: UM totales del plan x 2Gy/dosis por fracción en Gy
    MU_Gy = np.sum(MU_beam) * 2 / np.sum(d_beam) 
    return MU_Gy


def sas(data, umbral):
    umbral=float(input('Introducir umbral en mm'))#umbral por debajo del cual contamos un par de láminas como "poco abiertas".EN MILIMETROS!!!!

    acumulador_tot=0 #esto me sirve para contar el numero de pares de laminas abiertas en total (contando todos los haces)
    numero_cp_tot=0 #esto lo mismo pero con el numero de cp

    for beam in data.BeamSequence:
    
        for device in beam.BeamLimitingDeviceSequence:
            if device.RTBeamLimitingDeviceType=='MLCX':
                pares_laminas=device.NumberOfLeafJawPairs #cuento el número de pares de láminas 
        numero_cp=beam.NumberOfControlPoints #también el número de puntos de control (puede no ser el mismo para cada beam)
        print(beam.BeamName)
        print('Número de puntos de control: ',numero_cp)
        print('Número de pares de láminas:', pares_laminas)
    
        acumulador=0 #aqui contaré el numero de pares de laminas abiertas para cada haz singularmente
        for cp in beam.ControlPointSequence: #calculamos para cada punto de control (apertura) en cada beam
            
            for colimador in cp.BeamLimitingDevicePositionSequence: 
                
                if colimador.RTBeamLimitingDeviceType=='MLCX': 
                    posiciones=colimador.LeafJawPositions #asigno al vector posiciones las posiciones de todas las láminas
                    for i in range(pares_laminas):
                        d=posiciones[i+pares_laminas]-posiciones[i] #calculo la distancia a la que se encuentra cada par de láminas
                        if d > 0 and d < umbral: #si d=0 están cerradas y no me interesan para el cálculo
                            acumulador=acumulador+1 #numero de pares de láminas abierta más que el umbral
                    
        acumulador_tot=acumulador_tot+acumulador #sumo al acumulador total las correspondientes a cada beam
        numero_cp_tot=numero_cp_tot+numero_cp #idem con los puntos de control
        
        
        SAS=acumulador/pares_laminas/numero_cp*100 #calculo el SAS de cada beam como la proporcion de pares abiertas entre las totales(pares de laminas en el MLC * numero de puntos de control)
        print('SAS(',umbral,'mm) =',SAS,'%')
        

    SAS_tot=acumulador_tot/pares_laminas/numero_cp_tot*100
    print('PLAN COMPLETO')
    print('SAS(',umbral,'mm) =',SAS_tot,'%')

def pi(data):
    MU_todos_beams=[] #vector que almacena en el elemento i-ésimo las MU del beam i+1

    for fraction in data.FractionGroupSequence:
        for reference in fraction.ReferencedBeamSequence:
            MU_todos_beams.append(reference.BeamMeterset) #asigno a cada elemento las MU de cada beam por orden

            
    BI_beam=[] #aqui almacenaré el beam irregularity de cada beam
    PI_i=[] #esto es seria cada elemento del numerador del plan irregularity
    beam_number=0
    for beam in data.BeamSequence:
        beam_number=beam_number+1 #asigno un número de beam a cada uno.
        MU_beam=MU_todos_beams[beam_number-1] #digo cuantas MU tiene el beam con el que estamos trabajando
        
        FinalCumulativeMetersetWeight=beam.FinalCumulativeMetersetWeight #creo que siempre es 1, pero por si acaso
    
        for device in beam.BeamLimitingDeviceSequence:
            if device.RTBeamLimitingDeviceType=='MLCX':
                pares_laminas=device.NumberOfLeafJawPairs #cuento el número de pares de láminas 
                fronteras_MLC=device.LeafPositionBoundaries #las posiciones de los extremos de cada lámina.pares_laminas + 1 valores
        anchura=[]
        for i in range(pares_laminas):
            anchura.append(abs(fronteras_MLC[i]-fronteras_MLC[i+1])) #me calculo la anchura de cada lámina
            
        numero_cp=int(beam.NumberOfControlPoints) #también el número de puntos de control (puede no ser el mismo para cada beam)
        
        print(beam.BeamName,beam_number,MU_beam) #solo comprobacion para ver que todo marcha bien
        
        
        MU_cp_cumulative=[] #aquí iré metiendo el meterset cumulative para cada cp
        AI_cp=[] #aperture irregularity de cada cp
        
        for cp in beam.ControlPointSequence: #calculamos para cada punto de control (apertura) en cada beam
            MU_cp_cumulative.append(cp.CumulativeMetersetWeight*MU_beam/FinalCumulativeMetersetWeight)#cumulative meterset ya en UM para cada cp
            
            for colimador in cp.BeamLimitingDevicePositionSequence: 
                if colimador.RTBeamLimitingDeviceType=='MLCX': 
                    posiciones_izq=colimador.LeafJawPositions[:pares_laminas] #array con las posiciones de la parte izq del MLC
                    posiciones_der=colimador.LeafJawPositions[pares_laminas:] #array con las posiciones de la parte der del MLC
            
            # para el calculo del permitro de cada abertura hago dos pasos: primero asigno a cada par de laminas un índice m que es igual a 0
            # si ellas no pertenecen a ningún agujero. Por el contrario les asigno el número del agujero al que pertenecen.
            
            
            
            n_hole=0 # numero de agujeros en cada cp
            m=np.zeros(pares_laminas) # el elemento i-ésimo de este vector nos da m_i, que vale 0 si el par de láminas i-ésimo está cerrado ó vale el número del agujero al que dicho par de láminas pertence

            if posiciones_izq[0]==posiciones_der[0]: #primera lámina cerrada (hay que tomar el primer y último par de láminas por separado)
                m[0]=0
            else: #primera lámina abierta
                n_hole=1
                m[0]=n_hole

            for i in range(1,pares_laminas): #aqui corro desde el segundo par hasta el ultimo

                if posiciones_izq[i]==posiciones_der[i]: #lamina i-ésima cerrada
                    m[i]=0
                elif posiciones_izq[i]!=posiciones_der[i] and posiciones_izq[i-1]==posiciones_der[i-1]: #lamina i-ésima abierta y la anterior cerrada (empieza un nuevo agujero)
                    n_hole=n_hole+1
                    m[i]=n_hole
                elif posiciones_izq[i]!=posiciones_der[i] and posiciones_izq[i]>posiciones_der[i-1]: # la lamina actual (aun estando abierta) cierra el anterior agujero
                    n_hole=n_hole+1
                    m[i]=n_hole
                elif posiciones_izq[i]!=posiciones_der[i] and posiciones_der[i]<posiciones_izq[i-1]: #mismo caso que el anterior 
                    n_hole=n_hole+1
                    m[i]=n_hole
                else: #seguimos en el mismo agujero
                    m[i]=n_hole

            #Llegados a este punto ya tenemos, para cada cp, el índice de cada par de láminas.
            #vamos con el perimetro
            
            perimetro_cp=0 #perimetro total del cp (suma de los perimetros de cada agujero)
            perimetro_m=np.zeros(n_hole) #array donde cada elemento es el periemtro de cada uno de los agujeros en el cp
            
            if m[0]==0: #perimera lamina cerrada (m=0)
                pass
            elif m[1]!=m[0]: #la primera lamina es distinta que la segunda. Empieza y acaba un aguejro
                perimetro_m[0]=perimetro_m[0]+2*(posiciones_der[0]-posiciones_izq[0])+2*(anchura[0])

            else: #la primera es igual a la segunda. Empieza un agujero y sigue
                perimetro_m[0]=perimetro_m[0]+(posiciones_der[0]-posiciones_izq[0])+2*(anchura[0])

            for i in range(1,pares_laminas-1): #calculamos desde la segunda a la penultima

                if m[i]==0: #primera lamina cerrada
                    pass
                elif m[i-1]!=m[i] and m[i+1]==m[i]: #empieza un nuevo agujero y sigue
                    perimetro_m[int(m[i]-1)]=perimetro_m[int(m[i]-1)]+(posiciones_der[i]-posiciones_izq[i])+2*(anchura[i])
                elif m[i-1]!=m[i] and m[i+1]!=m[i]: #empieza un nuevo agujero y termina ahí
                    perimetro_m[int(m[i]-1)]=perimetro_m[int(m[i]-1)]+2*(posiciones_der[i]-posiciones_izq[i])+2*(anchura[i])
                elif m[i+1]==m[i]: #no empieza agujera (seguimos en el mismo que el anterior) y tampoco acaba
                    perimetro_m[int(m[i]-1)]=perimetro_m[int(m[i]-1)]+abs(posiciones_izq[i]-posiciones_izq[i-1])+abs(posiciones_der[i]-posiciones_der[i-1])+2*(anchura[i])
                else: #no empieza un agujero, pero si termina ahí
                    perimetro_m[int(m[i]-1)]=perimetro_m[int(m[i]-1)]+abs(posiciones_izq[i]-posiciones_izq[i-1])+abs(posiciones_der[i]-posiciones_der[i-1])+(posiciones_der[i]-posiciones_izq[i])+2*(anchura[i])
            
            # hago lo correspondiente para el último par de láminas
            if m[-1]==0: #ultima lamina cerrada
                perimetro_m[-1]=perimetro_m[-1]+(posiciones_der[-2]-posiciones_izq[-2])
            else: #ultima lamina abierta
                perimetro_m[-1]=perimetro_m[-1]+abs(posiciones_izq[-2]-posiciones_izq[-1])+abs(posiciones_der[-2]-posiciones_der[-1])+(posiciones_der[-1]-posiciones_izq[-1])+2*(anchura[-1])

            perimetro_cp=sum(perimetro_m) #perimetro de cada cp del beam en cuestion
            
            
            A_i=[] # lista donde el elemento i-ésimo es el area correspondiente al par de láminas i.
            for i in range(pares_laminas):
                A_i.append(anchura[i]*(posiciones_der[i]-posiciones_izq[i]))
            A_cp=sum(A_i)
            
            
            AI_cp.append(perimetro_cp**2/(4*np.pi*A_cp)) # aperture irregularity del cp
            
            

        #esta parte se dedica al calculo de las UM que le asigno a cada cp
        
        MU_cp=[] #vector vacio donde meteremos las MU que le asignamos a cada cp   
        MU_cp.append((MU_cp_cumulative[1]-MU_cp_cumulative[0])*0.5) #metemos la primera a mano. 
        for i in range(1,len(MU_cp_cumulative)-1):
            MU_cp.append((MU_cp_cumulative[i+1]-MU_cp_cumulative[i-1])*0.5) #expresión para las UM que le "asigno" a cada cp o apertura
        MU_cp.append((MU_cp_cumulative[-1]-MU_cp_cumulative[-2])*0.5) #la última igualmente a mano.
        
        MUij_AIij=[] #esto seria cada elemento del sumatorio del numerador del beam irregularity
        for i in range(numero_cp):
            MUij_AIij.append(AI_cp[i]*MU_cp[i])
        
        BI=sum(MUij_AIij)/MU_beam
        print('Beam Irregularity:',BI)
        
        
        
        PI_i.append(BI*MU_beam)
    PI=sum(PI_i)/sum(MU_todos_beams)

    print('Plan Completo')
    print('Plan Irregularity:',PI)