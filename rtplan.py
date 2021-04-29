import numpy as np
from pydicom.dataset import Dataset
from typing import Tuple, Dict

def getBeamLimitingDevice(name:str, node:Dataset) -> Dataset:
    """
    Busca un colimador de tipo 'name' en el nodo DICOM especificado 
    """
    lista = [x for x in node.BeamLimitingDeviceSequence if x.RTBeamLimitingDeviceType == name]
    if len(lista) > 1:
        raise ValueError(f"Se encontró más de un colimador {name} en el nodo")
    elif len(lista) == 0:
        raise ValueError(f"Colimador {name} no encontrado")
    return lista[0]

def getBeamLimitingDevicePosition(name:str, node:Dataset) -> Dataset:
    """
    Busca un posicionamiento de colimador de tipo 'name' en el nodo DICOM especificado 
    """
    lista = [x for x in node.BeamLimitingDevicePositionSequence if x.RTBeamLimitingDeviceType == name]
    if len(lista) > 1:
        raise ValueError(f"Se encontró más de un colimador {name} en el nodo")
    elif len(lista) == 0:
        raise ValueError(f"Colimador {name} no encontrado")
    return lista[0]

def get_mlc_geometry(nodo) -> Tuple[int, np.ndarray, np.ndarray]:
    """
    Devuelve la geometría del MLC definido en un nodo DICOM

    Parámetros
    ==========
    * nodo: DicomDataSet 
        El nodo con la información del MLC

    Returns
    =======
    * numero de pares de láminas
    * bordes de las láminas
    * anchuras de las láminas
    """
    mlc = getBeamLimitingDevice("MLCX", nodo)
    boundaries = np.array(mlc.LeafPositionBoundaries)
    widths = boundaries[1:] - boundaries[:-1]
    return mlc.NumberOfLeafJawPairs, boundaries, widths

def get_mlc_positions(nodo:Dataset) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Devuelve las posiciones del MLC definido en un nodo DICOM

    Parámetros
    ==========
    * nodo: DicomDataSet
        El nodo con la información del MLC

    Returns
    =======
    * posiciones de todas las láminas: array
    * posiciones de las láminas de la bancada izquierda: array
    * posiciones de las láminas de la bancada derecha: array
    """
    mlc = getBeamLimitingDevicePosition("MLCX", nodo)
    mlc_pos = np.array(mlc.LeafJawPositions)
    pares_laminas = len(mlc_pos) // 2
    posiciones_izq = mlc_pos[:pares_laminas] #array con las posiciones de la parte izq del MLC
    posiciones_der = mlc_pos[pares_laminas:] #array con las posiciones de la parte der del MLC
    return mlc_pos, posiciones_izq, posiciones_der


def get_beam_mu(data:Dataset) -> Dict[int, float]:
    """
    Devuelve un diccionario con el número de MU de cada haz
    """
    beam_mu = {}
    for fraction in data.FractionGroupSequence:
        for reference in fraction.ReferencedBeamSequence:
            beam_mu[reference.ReferencedBeamNumber] = reference.BeamMeterset
    return beam_mu