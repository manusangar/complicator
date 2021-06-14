
from typing import Tuple, List
import collections
import logging

import numpy as np
from pydicom.dataset import Dataset

from rtplan import RTPlan
from rtplan import get_mlc_geometry, get_mlc_positions, get_beam_mu

def mu_per_gy(plan:RTPlan) -> float:
    """
    Devuelve el índice de complejidad MU/Gy para el plan especificado en data
    """
    total_mu = sum([beam.mu for beam in plan.beams.values()])
    total_dose = sum([beam.dose for beam in plan.beams.values()])
    MU_Gy = total_mu  * 2 / total_dose
    return MU_Gy

def sas(plan:RTPlan, umbral:float) -> float:
    """
    Calcula el índice SAS de complejidad

    El SAS se define como la proporcion de pares de láminas abiertas menos
    de un cierto umbral con respecto al número de láminas total

    Parámetros
    ==========
    * plan: RTPlan
        El plan a evaluar
    * umbral:float
        El umbral por debajo del cual consideramos que un par de láminas están poco abiertas (en mm)
    """
    total_open_leaves = 0 #esto me sirve para contar el numero de pares de laminas abiertas en total (contando todos los haces)
    numero_cp_tot = 0 #esto lo mismo pero con el numero de cp

    for beam_id, beam in plan.beams.items():

        numero_cp_tot += beam.number_of_segments #también el número de puntos de control (puede no ser el mismo para cada beam)
        beam_open_leaves = 0
        logging.debug("SAS beam %d, beam, cp, open_leaves", beam_id) 
        for idx, segment in enumerate(beam.segments): #calculamos para cada punto de control (apertura) en cada beam
            d = segment.mlc_right - segment.mlc_left
            cp_open_leaves = np.count_nonzero(np.logical_and(d > 0, d < umbral))
            beam_open_leaves += cp_open_leaves
            logging.debug(f"SAS cp, %d, %d, %d",beam_id, idx, cp_open_leaves)                    
        total_open_leaves += beam_open_leaves 
 
        logging.debug(f"SAS beam: %d %.2f", beam_id, 
                                            beam_open_leaves / (beam.mlc_geometry.leaf_number * beam.number_of_segments) * 100)

    SAS_tot = 100 * total_open_leaves / (beam.mlc_geometry.leaf_number * numero_cp_tot)
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


def pi(plan:RTPlan) -> float:
    """
    Calcula el índice de complejidad PI para un plan de VMAT/IMR

    Parámetros
    ==========
    * data: DicomDataset
        El plan de tratamiento en formato DICOM
    """            
    PI_i = [] #esto es seria cada elemento del numerador del plan irregularity

    total_mu = 0
    for beam_id, beam in plan.beams.items():
        total_mu += beam.mu
        logging.debug("PI beam start %d", beam_id)
       
        AI_cp = np.zeros(beam.number_of_segments) #aperture irregularity de cada cp
        BI_cp = np.zeros(beam.number_of_segments) #aperture irregularity de cada cp
        logging.debug("PI cp, beam, cp, perimetro, apertura, irregularidad")
        for cp_idx, segment in enumerate(beam.segments): #calculamos para cada punto de control (apertura) en cada beam
            perimetro_cp = get_perimetro(segment.mlc_left, segment.mlc_right, beam.mlc_geometry.leaf_widths)
            A_cp = np.sum(beam.mlc_geometry.leaf_widths * (segment.mlc_right - segment.mlc_left))
            AI_cp[cp_idx] = perimetro_cp**2 / (4 * np.pi * A_cp) # aperture irregularity del cp
            BI_cp[cp_idx] = (AI_cp[cp_idx] * segment.mu)
            logging.debug("PI cp, %d, %d, %.2f, %.2f, %.2f", beam_id, cp_idx, perimetro_cp, A_cp, AI_cp[cp_idx])
        
        BI = np.sum(BI_cp) / beam.mu
        PI_i.append(BI * beam.mu)
        logging.debug("PI beam, %d, MU = %.2f, BI = %.2f", beam_id, beam.mu, BI)

    PI= sum(PI_i) / total_mu
    return PI


def tgi(plan:RTPlan) -> float:
    "Calcula el indice Tongue and Groove"

    total_mu = sum([beam.mu for beam in plan.beams.values()])
    mean_gap = 0
    mean_tgi = 0
    for beam_id, beam in plan.beams.items():
        for segment in beam.segments:
            
            gap = segment.mlc_right - segment.mlc_left
            overlap = np.abs(segment.mlc_left[1:] - segment.mlc_left[:-1]) + np.abs(segment.mlc_right[1:] - segment.mlc_right[:-1])

            tgi = np.zeros(len(overlap))
            open_leafs = np.logical_not(np.isclose(segment.mlc_right, segment.mlc_left))[1:]
            tgi_gap = gap[1:]
            tgi[open_leafs] = np.minimum(overlap[open_leafs] / tgi_gap[open_leafs], 1.0)
            mean_tgi += np.mean(tgi) * segment.mu 
            mean_gap += np.mean(gap) * segment.mu
            logging.debug("TGI %f, %f", mean_tgi, mean_gap)
    mean_gap /= total_mu
    mean_tgi /= total_mu
    return mean_gap, mean_tgi
