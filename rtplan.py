import numpy as np

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

def get_mlc_geometry(nodo):
    """
    Devuelve la geometría del MLC definido en un nodo DICOM

    Parameters
    ==========
    * nodo: DicomDataSet. El nodo con la información del MLC

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

def get_mlc_positions(nodo):
    """
    Devuelve las posiciones del MLC definido en un nodo DICOM

    Parameters
    ==========
    * nodo: DicomDataSet. El nodo con la información del MLC

    Returns
    =======
    * posiciones de todas las láminas
    * posiciones de las láminas de la bancada izquierda
    * posiciones de las láminas de la bancada derecha
    """
    mlc = getBeamLimitingDevicePosition("MLCX", nodo)
    mlc_pos = np.array(mlc.LeafJawPositions)
    pares_laminas = len(mlc_pos) // 2
    posiciones_izq = mlc_pos[:pares_laminas] #array con las posiciones de la parte izq del MLC
    posiciones_der = mlc_pos[pares_laminas:] #array con las posiciones de la parte der del MLC
    return mlc_pos, posiciones_izq, posiciones_der


def get_beam_mu(data):
    """
    Devuelve un diccionario con el número de MU de cada haz
    """
    beam_mu = {}
    for fraction in data.FractionGroupSequence:
        for reference in fraction.ReferencedBeamSequence:
            beam_mu[reference.ReferencedBeamNumber] = reference.BeamMeterset
    return beam_mu