
from typing import Tuple, List
import collections
import logging

import numpy as np
from pydicom.dataset import Dataset

from rtplan import get_mlc_geometry, get_mlc_positions, get_beam_mu

def mu_per_gy(data:Dataset) -> float:
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

def sas(data:Dataset, umbral:float) -> float:
    """
    Calcula el índice SAS de complejidad

    El SAS se define como la proporcion de pares de láminas abiertas menos
    de un cierto umbral con respecto al número de láminas total

    Parámetros
    ==========
    * data: DicomDataset
        El objeto DICOM con el plan a evaluar
    * umbral:float
        El umbral por debajo del cual consideramos que un par de láminas están poco abiertas (en mm)
    """
    total_open_leaves = 0 #esto me sirve para contar el numero de pares de laminas abiertas en total (contando todos los haces)
    numero_cp_tot = 0 #esto lo mismo pero con el numero de cp

    for beam in data.BeamSequence:
        pares_laminas, _, _  = get_mlc_geometry(beam)

        numero_cp_tot += beam.NumberOfControlPoints #también el número de puntos de control (puede no ser el mismo para cada beam)
        beam_open_leaves = 0 
        for cp in beam.ControlPointSequence: #calculamos para cada punto de control (apertura) en cada beam
            _, posiciones_izq, posiciones_dcha = get_mlc_positions(cp)
            d = posiciones_dcha - posiciones_izq
            cp_open_leaves = np.count_nonzero(np.logical_and(d > 0, d < umbral))
            beam_open_leaves += cp_open_leaves
            #print(f"DEBUG SAS cp: {beam.BeamNumber}, {cp.ControlPointIndex}, {cp_open_leaves}")                    
        total_open_leaves += beam_open_leaves 

        #SAS = beam_open_leaves / (pares_laminas * numero_cp) * 100 
        #print(f"DEBUG SAS beam: {beam.BeamNumber}, {SAS:.2f}%")

    SAS_tot = 100 * total_open_leaves / (pares_laminas * numero_cp_tot)
    return SAS_tot


class MLCHole:
    """
    Describe un agujero en el MLC identificado como un par de
    número de láminas, correspondientes a la primera y última
    abiertas
    """
    def __init__(self, first_leaf, last_leaf):
        self.first_leaf = first_leaf
        self.last_leaf = last_leaf


def find_holes(posiciones_izq:np.ndarray, posiciones_der:np.ndarray) -> List[MLCHole]:
    """
    Busca los agujeros en una conformación del MLC

    Parámetros
    ==========
    * posiciones_izq: array
        Array con las posiciones de las láminas de la bancada izquierda
    * posiciones_dcha: array
        Array con las posiciones de las láminas de la bancada derecha

    Returns
    =======
    Array de objetos MLCHole con los agujeros del MLC
    """
    holes = []
    current_hole = None
    leaf_closed = np.isclose(posiciones_izq, posiciones_der)

    for i in range(len(posiciones_izq)): #aqui corro desde el segundo par hasta el ultimo
        if leaf_closed[i]:
            #Agujero cerrado
            current_hole = None
        else:
            if current_hole is None or posiciones_izq[i] > posiciones_der[i-1] or posiciones_der[i] < posiciones_izq[i-1]:
                #Empieza un nuevo agujero porque el anterior estaba cerrado o porque
                #las láminas cierran el anterior a la vez que abren el nuevo
                current_hole = MLCHole(i,i)
                holes.append(current_hole)
            else:
                current_hole.last_leaf = i
    return holes


def get_perimetro(posiciones_izq:np.ndarray, posiciones_der:np.ndarray, anchura:np.ndarray) -> float:
    """
    Calcula el perímetro de las aperturas en una conformación del MLC

    Parámetros
    ==========
    * posiciones_izq: array
        Array con las posiciones de las láminas de la bancada izquierda
    * posiciones_dcha: array
        Array con las posiciones de las láminas de la bancada derecha
    * anchura: array
        Array con las anchuras de las láminas del MLC

    Returns
    =======
    Suma de los perímetros de los agujeros del MLC
    """
    # para el calculo del permitro de cada abertura hago dos pasos: primero asigno a cada par de laminas un índice m que es igual a 0
    # si ellas no pertenecen a ningún agujero. Por el contrario les asigno el número del agujero al que pertenecen.  
    pares_laminas = len(posiciones_izq)

    holes = find_holes(posiciones_izq, posiciones_der)

    perimetro_cp = 0
    for hole in holes:
        #Borde superior
        hole_perimeter = posiciones_der[hole.first_leaf] - posiciones_izq[hole.first_leaf]

        for leaf in range(hole.first_leaf + 1, hole.last_leaf + 1):
            hole_perimeter += abs(posiciones_izq[leaf] - posiciones_izq[leaf - 1]) + \
                              abs(posiciones_der[leaf] - posiciones_der[leaf - 1])

        #Borde inferior
        hole_perimeter += posiciones_der[hole.last_leaf] - posiciones_izq[hole.last_leaf]
        hole_perimeter += 2 * np.sum(anchura[hole.first_leaf:hole.last_leaf + 1])

        perimetro_cp += hole_perimeter

    return perimetro_cp


def pi(data:Dataset) -> float:
    """
    Calcula el índice de complejidad PI para un plan de VMAT/IMR

    Parámetros
    ==========
    * data: DicomDataset
        El plan de tratamiento en formato DICOM
    """
    beam_mu = get_beam_mu(data)
            
    PI_i = [] #esto es seria cada elemento del numerador del plan irregularity

    for beam in data.BeamSequence:
        logging.debug("PI beam start %d", beam.BeamNumber)
        MU_beam = beam_mu[beam.BeamNumber] #digo cuantas MU tiene el beam con el que estamos trabajando
        _, _, anchuras = get_mlc_geometry(beam)
            
        numero_cp=int(beam.NumberOfControlPoints) #también el número de puntos de control (puede no ser el mismo para cada beam)
        
        print(beam.BeamName, beam.BeamNumber ,beam_mu[beam.BeamNumber]) #solo comprobacion para ver que todo marcha bien
            
        cp_cumulative_weight = np.zeros(numero_cp)
        AI_cp = np.zeros(numero_cp) #aperture irregularity de cada cp
        logging.debug("PI cp, beam, cp, perimetro, apertura, irregularidad")
        for cp_idx, cp in enumerate(beam.ControlPointSequence): #calculamos para cada punto de control (apertura) en cada beam
            cp_cumulative_weight[cp_idx] = cp.CumulativeMetersetWeight
            
            _, posiciones_izq, posiciones_der = get_mlc_positions(cp)
            
            perimetro_cp = get_perimetro(posiciones_izq, posiciones_der, anchuras)
            A_cp = np.sum(anchuras * (posiciones_der - posiciones_izq))
            AI_cp[cp_idx] = perimetro_cp**2 / (4 * np.pi * A_cp) # aperture irregularity del cp
            logging.debug("PI cp, %d, %d, %.2f, %.2f, %.2f", beam.BeamNumber, cp_idx, perimetro_cp, A_cp, AI_cp[cp_idx])
            
        #esta parte se dedica al calculo de las UM que le asigno a cada cp
        MU_cp_cumulative = MU_beam * cp_cumulative_weight / beam.FinalCumulativeMetersetWeight
        MU_cp = (np.diff(MU_cp_cumulative, append=MU_beam) + np.diff(MU_cp_cumulative, prepend=0)) / 2.0
        
        BI= np.sum(AI_cp * np.array(MU_cp)) / MU_beam
        PI_i.append(BI * MU_beam)
        logging.debug("PI beam, %d, MU = %.2f, BI = %.2f", beam.BeamNumber, MU_beam, BI)

    PI= sum(PI_i) / sum(beam_mu.values())
    return PI