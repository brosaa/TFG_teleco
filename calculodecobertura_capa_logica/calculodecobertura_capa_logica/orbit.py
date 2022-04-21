"""calculodecobertura_capa_logica - Capa lógica o de negocios del software de cálculo de cobertura de un satélite, Trabajo de Fin de Grado en ingeniería de tecnologías de las telecomunicaciones URJC 2022 de Belén Rosa Alonso"""

__version__ = '0.0.1'
__author__ = 'Belén Rosa Alonso <b.rosaa.2017@alumnos.urjc.es>'
__all__ = []

from math import pi, sin, asin
from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import propagate
from poliastro.bodies import Earth

from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)

def visionArc(Rtierra, alt, FOV):

    lado1 = Rtierra
    lado2 = (alt / u.km) + Rtierra

    angulo1 = FOV/2
    angulo1_rad = (angulo1* pi)/180

    angulo2 = 180 - (asin(lado2 * sin(angulo1_rad) / lado1) * 180/pi)

    angulo3_rad = (180 - angulo1 - angulo2) * pi/180

    return Rtierra * angulo3_rad

def createOrbit(Earth, alt, inc, raan):
    sat = Orbit.circular(Earth, alt, inc, raan)

    return sat

def orbitTimes(cycle, timeDelta):
    """
    Función que devuelve una lista con los tiempos desde un instante inicial, hasta un rango de minutos después
    """
    orbitRange = cycle * 60
    time_init = "00:00"

    orbitTimes = []

    hour = time_init.split(':')[0]
    minute = time_init.split(':')[1]

    orbitTimes.append(time_init)
 
    for i in range(0, orbitRange, timeDelta):
        if timeDelta < 60:
            minute = int(minute) + timeDelta

            if minute >= 60:
                hour = int(hour) + 1                
                minute = minute - 60 
        elif timeDelta == 60:
            hour = int(hour) + 1
            minute = int(minute)
        else:
            hours_plus = timeDelta / 60
            hour = int(hour) + int(hours_plus)
            minute = int(minute) + int((hours_plus - int(hours_plus)) * 60)

            if minute >= 60:
                minute = minute - 60
                hour = hour + 1
        
        if minute <= 9:
            minute = "0" + str(minute)
        
        if int(hour) <= 9:
            hour = "0" + str(int(hour))
        elif int(hour) == 24:
            hour = "00"

        time = str(hour) + ":" + str(minute)
        orbitTimes.append(time)
    
    return orbitTimes

def _get_raw_coords(orb, t_deltas):
    """Generates raw orbit coordinates for given epochs

    Parameters
    ----------
    orb: ~poliastro.twobody.orbit
        Orbit to be propagated
    t_deltas: ~astropy.time.DeltaTime
        Desired observation time

    Returns
    -------
    raw_xyz: array
        A collection of raw cartessian position vectors
    raw_epochs: array
        Associated epoch with previously raw coordinates
    """

    # Solve for raw coordinates and epochs
    raw_xyz = propagate(orb, t_deltas)
    raw_epochs = orb.epoch + t_deltas

    return raw_xyz, raw_epochs

def _from_raw_to_ITRS(raw_xyz, raw_obstime):
    """Converts raw coordinates to ITRS ones

    Parameters
    ----------
    raw_xyz: array
        A collection of rwa position coordinates
    raw_obstime: array
        Associated observation time

    Returns
    -------
    itrs_xyz: ~astropy.coordinates.ITRS
        A collection of coordinates in ITRS frame

    """

    # Build GCRS and ITRS coordinates
    gcrs_xyz = GCRS(
        raw_xyz, obstime=raw_obstime, representation_type=CartesianRepresentation
    )
    itrs_xyz = gcrs_xyz.transform_to(ITRS(obstime=raw_obstime))

    return itrs_xyz       

def sphCoords(orb, t_deltas):
    """
    Esta función convierte coordenadas cartesianas (x,y,z) en coordenadas
    esféricas, formato necesario para otros programas que usan este como
    entrada. Hay que tener en cuenta que el plano ecuatorial está formado
    por los ejes "x,y", mientras que el eje z está contenido en el eje de 
    rotación de la tierra
    """
    # Compute predicted grountrack positions
    raw_xyz, raw_obstime = _get_raw_coords(orb, t_deltas)
    itrs_xyz = _from_raw_to_ITRS(raw_xyz, raw_obstime)
    itrs_latlon = itrs_xyz.represent_as(SphericalRepresentation)

    # Append predicted positions to map
    lat=itrs_latlon.lat.to(u.deg)/ u.deg
    lon=itrs_latlon.lon.to(u.deg)/ u.deg
    return [float(lat), float(lon) - 180]

def propagateSat(times, timeDelta, alt, inc, raan):
    polar_positions = {}
    propagationTimeDelta_init = timeDelta * u.min
    propagationTimeDelta = 0 * u.min
    
    #Creamos la órbita del satélite final con hora correcta
    orb = createOrbit(Earth, alt, inc, raan)

    print("Calculating satellite positions...")

    for time in times:
        polar_positions[time] = sphCoords(orb, propagationTimeDelta)
        propagationTimeDelta = propagationTimeDelta + propagationTimeDelta_init
    
    return polar_positions