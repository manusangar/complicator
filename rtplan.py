import numpy as np
from pydicom.dataset import Dataset
from typing import Tuple, Dict

class RTPlan:
    def __init__(self, beams):
        self.beams = beams

    @classmethod
    def load_dicom(cls, data: Dataset):
        beam_mu, beam_dose = get_beam_mu(data)
        beams = {}
        for beam in data.BeamSequence:
            numero_cp = beam.NumberOfControlPoints
            MU_beam = beam_mu[beam.BeamNumber]
            Dose_beam = beam_dose[beam.BeamNumber]
            leaf_number, boundaries, widths = get_mlc_geometry(beam)
            segments = []

            cp_cumulative_weight = np.zeros(numero_cp)
            segments = []
            for cp_idx, cp in enumerate(beam.ControlPointSequence):
                cp_cumulative_weight[cp_idx] = cp.CumulativeMetersetWeight
                mlc_pos, posiciones_izq, posiciones_dcha = get_mlc_positions(cp)
                segments.append(RTSegment(posiciones_izq, posiciones_dcha, 0.0))
                
            MU_cp_cumulative = MU_beam * cp_cumulative_weight / beam.FinalCumulativeMetersetWeight
            MU_cp = (np.diff(MU_cp_cumulative, append=MU_beam) + np.diff(MU_cp_cumulative, prepend=0)) / 2.0
            for idx in range(len(MU_cp)):
                segments[idx].mu = MU_cp[idx]

            beams[beam.BeamNumber] = RTBeam(segments, MU_beam, Dose_beam, MLCGeometry( leaf_number, boundaries, widths))
        return cls(beams)

class MLCGeometry:
    def __init__(self, number_of_leafs, leaf_boundaries, leaf_widths):
        self.leaf_number = number_of_leafs
        self.leaf_boundaries = leaf_boundaries
        self.leaf_widths = leaf_widths

class RTBeam:
    def __init__(self, segments, mu, dose, mlc_geometry):
        self.segments = segments
        self.number_of_segments = len(segments)
        self.mu = mu
        self.dose = dose
        self.mlc_geometry = mlc_geometry

class RTSegment:
    def __init__(self, mlc_left, mlc_right, mu):
        self.mlc_left = mlc_left
        self.mlc_right = mlc_right
        self.mu = mu

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
    beam_dose = {}
    for fraction in data.FractionGroupSequence:
        for reference in fraction.ReferencedBeamSequence:
            beam_mu[reference.ReferencedBeamNumber] = reference.BeamMeterset
            beam_dose[reference.ReferencedBeamNumber] = reference.BeamDose
    return beam_mu, beam_dose