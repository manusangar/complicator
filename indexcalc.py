
from typing import NamedTuple

import numpy as np

def mu_per_gy(data):
    """
    Devuelve el índice de complejidad MU/Gy para el plan especificado en data
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

def getBeamLimitingDevice(name:str, node):
    lista = [x for x in node.BeamLimitingDeviceSequence if x.RTBeamLimitingDeviceType == name]
    if len(lista) > 1:
        raise ValueError(f"Se encontró más de un colimador {name} en el nodo")
    elif len(lista) == 0:
        raise ValueError(f"Colimador {name} no encontrado")
    return lista[0]

def getBeamLimitingDevicePosition(name:str, node):
    lista = [x for x in node.BeamLimitingDevicePositionSequence if x.RTBeamLimitingDeviceType == name]
    if len(lista) > 1:
        raise ValueError(f"Se encontró más de un colimador {name} en el nodo")
    elif len(lista) == 0:
        raise ValueError(f"Colimador {name} no encontrado")
    return lista[0]

def sas(data, umbral):
    """
    Calcula el índice SAS de complejidad

    El SAS se define como la proporcion de pares de láminas abiertas menos
    de un cierto umbral con respecto al número de láminas total

    Parameters
    ==========
    * umbral:float
        El umbral por debajo del cual consideramos que un par de láminas están poco abiertas (en mm)
    """
    total_open_leaves = 0 #esto me sirve para contar el numero de pares de laminas abiertas en total (contando todos los haces)
    numero_cp_tot = 0 #esto lo mismo pero con el numero de cp

    for beam in data.BeamSequence:
        mlc_desc = getBeamLimitingDevice("MLCX", beam)
        pares_laminas = mlc_desc.NumberOfLeafJawPairs

        numero_cp = beam.NumberOfControlPoints #también el número de puntos de control (puede no ser el mismo para cada beam)
        beam_open_leaves = 0 
        for cp in beam.ControlPointSequence: #calculamos para cada punto de control (apertura) en cada beam
            mlc_pos = getBeamLimitingDevicePosition("MLCX", cp)
            posiciones = np.array(mlc_pos.LeafJawPositions)
            d = posiciones[pares_laminas:] - posiciones[:pares_laminas]
            cp_open_leaves = np.count_nonzero(np.logical_and(d > 0, d < umbral))
            beam_open_leaves += cp_open_leaves
            #print(f"DEBUG SAS cp: {beam.BeamNumber}, {cp.ControlPointIndex}, {cp_open_leaves}")                    
        total_open_leaves += beam_open_leaves 
        numero_cp_tot += numero_cp 

        #SAS = beam_open_leaves / (pares_laminas * numero_cp) * 100 
        #print(f"DEBUG SAS beam: {beam.BeamNumber}, {SAS:.2f}%")

    SAS_tot = 100 * total_open_leaves / (pares_laminas * numero_cp_tot)
    return SAS_tot


def get_mlc_geometry(nodo):
    mlc = getBeamLimitingDevice("MLCX", nodo)
    boundaries = np.array(mlc.LeafPositionBoundaries)
    widths = boundaries[1:] - boundaries[:-1]
    return mlc.NumberOfLeafJawPairs, boundaries, widths


def get_beam_mu(data):
    """
    Devuelve un diccionario con el número de MU de cada haz
    """
    beam_mu = {}
    for fraction in data.FractionGroupSequence:
        for reference in fraction.ReferencedBeamSequence:
            beam_mu[reference.ReferencedBeamNumber] = reference.BeamMeterset
    return beam_mu

def get_perimetro(posiciones_izq, posiciones_der, anchura):
    # para el calculo del permitro de cada abertura hago dos pasos: primero asigno a cada par de laminas un índice m que es igual a 0
    # si ellas no pertenecen a ningún agujero. Por el contrario les asigno el número del agujero al que pertenecen.  
    pares_laminas = len(posiciones_izq)
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
    return perimetro_cp


def pi(data):
    beam_mu = get_beam_mu(data)
            
    PI_i = [] #esto es seria cada elemento del numerador del plan irregularity

    for beam in data.BeamSequence:
        MU_beam = beam_mu[beam.BeamNumber] #digo cuantas MU tiene el beam con el que estamos trabajando
        pares_laminas, _, anchura = get_mlc_geometry(beam)
            
        numero_cp=int(beam.NumberOfControlPoints) #también el número de puntos de control (puede no ser el mismo para cada beam)
        
        print(beam.BeamName, beam.BeamNumber ,beam_mu[beam.BeamNumber]) #solo comprobacion para ver que todo marcha bien
        
        
        MU_cp_cumulative=[] #aquí iré metiendo el meterset cumulative para cada cp
        AI_cp=[] #aperture irregularity de cada cp
        
        for cp in beam.ControlPointSequence: #calculamos para cada punto de control (apertura) en cada beam
            MU_cp_cumulative.append(cp.CumulativeMetersetWeight*MU_beam/beam.FinalCumulativeMetersetWeight)#cumulative meterset ya en UM para cada cp
            
            mlc_cp = getBeamLimitingDevicePosition("MLCX", cp)
            mlc_cp_pos = np.array(mlc_cp.LeafJawPositions)
            posiciones_izq = mlc_cp_pos[:pares_laminas] #array con las posiciones de la parte izq del MLC
            posiciones_der = mlc_cp_pos[pares_laminas:] #array con las posiciones de la parte der del MLC
            
            perimetro_cp = get_perimetro(posiciones_izq, posiciones_der, anchura)

            A_cp = np.sum(anchura * (posiciones_der - posiciones_izq))
          
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
        
        BI= sum(MUij_AIij) / beam_mu[beam.BeamNumber]
        print('Beam Irregularity:',BI)
        
        
        
        PI_i.append(BI * beam_mu[beam.BeamNumber])
    PI= sum(PI_i) / sum(beam_mu.values())

    print('Plan Completo')
    print('Plan Irregularity:',PI)
    return PI