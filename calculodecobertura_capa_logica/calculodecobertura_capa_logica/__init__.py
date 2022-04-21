"""calculodecobertura_capa_logica - Capa lógica o de negocios del software de cálculo de cobertura de un satélite, Trabajo de Fin de Grado en ingeniería de tecnologías de las telecomunicaciones URJC 2022 de Belén Rosa Alonso"""

__version__ = '0.0.1'
__author__ = 'Belén Rosa Alonso <b.rosaa.2017@alumnos.urjc.es>'
__all__ = []

import math
from astropy import units as u

def visionArc(Rtierra, alt, FOV):

    lado1 = Rtierra
    lado2 = (alt / u.km) + Rtierra

    angulo1 = FOV/2
    angulo1_rad = (angulo1* math.pi)/180

    angulo2 = 180 - (math.asin(lado2 * math.sin(angulo1_rad) / lado1) * 180/math.pi)

    angulo3_rad = (180 - angulo1 - angulo2) * math.pi/180
    return Rtierra * angulo3_rad
