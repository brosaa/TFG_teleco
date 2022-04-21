
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