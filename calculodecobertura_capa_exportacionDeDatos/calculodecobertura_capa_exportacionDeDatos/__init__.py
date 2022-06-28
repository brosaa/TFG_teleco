"""calculodecobertura_capa_exportacionDeDatos - Capa de exportación de los datos resultado del software de cálculo de cobertura de un satélite, Trabajo de Fin de Grado en ingeniería de tecnologías de las telecomunicaciones URJC 2022 de Belén Rosa Alonso"""

__version__ = '0.0.1'
__author__ = 'Belén Rosa Alonso <b.rosaa.2017@alumnos.urjc.es>'
__all__ = []

from csv import writer

def exportObjectives(objectivesInFootprint, objectives):
    
    # Nombre del archivo final
    out_file = "./objectivesInFootprint.csv"
    print("Exporting objectives under the satellite footprint at " + out_file + " file...")

    data = []
    header = []

    for objective in objectives:
        header.append(objective)
    
    header.append("Visual precision (%)")
    header.append("Satellite position lat")
    header.append("Satellite position lon")
    header.append("Time since AN")

    data.append(header)


    for objective in objectivesInFootprint:

        objective_data = [objective]

        for i in range(len(objectivesInFootprint[objective]) - 2): 
            objective_data.append(objectivesInFootprint[objective][i])

        sat_pos_lat = objectivesInFootprint[objective][len(objectivesInFootprint[objective]) - 2][0]
        sat_pos_lon = objectivesInFootprint[objective][len(objectivesInFootprint[objective]) - 2][1]
        Time = objectivesInFootprint[objective][len(objectivesInFootprint[objective]) - 1]

        objective_data.append(sat_pos_lat)
        objective_data.append(sat_pos_lon)
        objective_data.append(Time)
        
        data.append(objective_data)

    csvsalida = open(out_file, 'w', newline='', encoding="utf-8")
    salida = writer(csvsalida)
    salida.writerows(data)
    del salida
    csvsalida.close()

def exportSatPositions(polar_positions):
    """
    Esta función almacena en un archivo Excel los datos referentes a si los FIR están
    cubiertos entera o parcialmente, o si no han sido cubiertos, al igual que el número
    de puntos cubiertos respecto al número de puntos inicial.
    """

    # Nombre del archivo final
    out_file = "./satellite_positions.csv"

    print("Exporting satellite positions to " + out_file + " file...")

    data = []
    header = ["Time since AN","lat","lon"]
    data.append(header)


    for time in polar_positions:

        time_data = [time] 
        lat = polar_positions[time][0]
        lon = polar_positions[time][1]

        time_data.append(lat)
        time_data.append(lon)
        
        data.append(time_data)

    csvsalida = open(out_file, 'w', newline='')
    salida = writer(csvsalida)
    salida.writerows(data)
    del salida
    csvsalida.close()