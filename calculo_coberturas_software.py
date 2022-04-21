#TFG Thales Alenia Space
# Orbital Mechanics

#Autoría: Belén Rosa Alonso
#Titulación: INGENIERÍA EN TECNOLOGÍAS DE LA TELECOMUNICACION + INGENIERÍA AEROESPACIAL EN AERONAVEGACION
#Tutores: Ángel Alvaro Sánchez y Rafael Jose Navarro Sebastiá
#TFG de cálculo de coberturas - Versión 2

#PAQUETES A IMPORTAR

#Paquetes propios
from calculodecobertura_capa_logica.prueba1 import visionArc
from calculodecobertura_capa_presentacion import close, checkValues, init

#Paquetes públicos
from astropy import units as u

from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import propagate
from astropy.coordinates import (
    GCRS,
    ITRS,
    CartesianRepresentation,
    SphericalRepresentation,
)

from tkinter import *
from shapely.geometry.polygon import Polygon, LinearRing, Point, LineString
import math
import great_circle_calculator.great_circle_calculator as gcc
import pandas as pd
import numpy as np
import csv
from sys import exit

#VARIABLES FIJAS

#Órbita circular
ecc = 0 * u.one

#Radio de la Tierra (km)
Rtierra = 6371

#Tiempo de propagación [min]
timeDelta = 3


#VARIABLES CONFIGURABLES

#Altitud orbital
alt = 562 * u.km

#Inclinación orbital
inc = 97.64 * u.deg

#Right ascension of the ascendint node
raan = 0 * u.deg

#Campo de visión del instrumento ACT
FOV = 10.19

#Ciclo de repetición del satélite en horas
cycle = 24


#Obtener coordenadas de las ciudades del mundo desde el excel
file = ""
objectives_init = ""

def checkObjetives(obj, data):
    if(not (obj in data)):
        print("Not valid objectives characteristics to calculate, please try again inserting valid arguments.\n")
        print("Make sure the characteristics are the column names from the input excel file.\n\nClosing program...")

        exit()

def getObjectives(excel, objectives):
    """
    Función para leer los datos obtenidos desde el excel de objetivos en Tierra
    """

    data_info = {}

    data = pd.read_excel(excel)

    for objective in objectives:
        if(objective[0] == " "):
            obj = objective[1:]   
        else:
            obj = objective

        checkObjetives(obj, data)
        data_info[objective] = data[obj].to_list()
        
    objs = {}
    
    for i in range(len(data_info[objectives[0]])):
        objs[data_info[objectives[0]][i]] = []
        for j in range(len(objectives) - 1):
            objs[data_info[objectives[0]][i]].append(data_info[objectives[j + 1]][i])

    return objs

def get_parameters():
    global alt, inc, raan, FOV, file, objectives_init, cycle

    if e1.get() != "":
        alt = int(e1.get()) * u.km

        actual_altitude = "Current altitude is " + str(alt) + ", introduce new altitude: "
        altitude.configure(text = actual_altitude)
    
    if e2.get() != "":
        inc = float(e2.get()) * u.deg

        actual_inclination = "Current inclination is " + str(inc) + ", introduce new inclination: "
        inclination.configure(text = actual_inclination)

    if e3.get() != "":
        raan = float(e3.get()) * u.deg

        actual_raan = "Current right ascension of the ascendint node is " + str(raan) + ", introduce new raan: "
        raan2.configure(text = actual_raan)

    if e4.get() != "":
        FOV = float(e4.get())

        actual_FOV = "Current FOV is " + str(FOV) + " degrees, introduce new FOV: "
        FOV2.configure(text = actual_FOV)

    if e5.get() != "":
        file = "./" + e5.get()

        if file == "":
            actual_fileinput = "Currently there's not file with objectives, introduce file name (example: worldcities.xlsx): "
        else:
            actual_fileinput = "Current Excel with objectives is " + str(file) + ", introduce new file name: "
        
        fileinput.configure(text = actual_fileinput)

    if e6.get() != "":
        objectives_init = e6.get().split(',')

        if objectives_init == "":
            actual_objectives = "Currently there're not objective's characteristics, introduce name, lat, lon and aditional characteristic: "
        else:
            actual_objectives = "Current objective's characteristics are " + str(objectives_init) + ", introduce new objective's characteristics: "

        objectives2.configure(text = actual_objectives)

    if e7.get() != "":
        cycle = int(e7.get())

        actual_cycle = "Current calculation period is " + str(cycle) + " hours, introduce new calculation period: "
        cycle2.configure(text = actual_cycle)

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

def propagateSat(times, orb):
    polar_positions = {}
    propagationTimeDelta_init = timeDelta * u.min
    propagationTimeDelta = 0 * u.min

    print("Calculating satellite positions...")

    for time in times:
        polar_positions[time] = sphCoords(orb, propagationTimeDelta)
        propagationTimeDelta = propagationTimeDelta + propagationTimeDelta_init
    
    return polar_positions
    
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
    header = ["Time","lat","lon"]
    data.append(header)


    for time in polar_positions:

        time_data = [time] 
        lat = polar_positions[time][0]
        lon = polar_positions[time][1]

        time_data.append(lat)
        time_data.append(lon)
        
        data.append(time_data)

    csvsalida = open(out_file, 'w', newline='')
    salida = csv.writer(csvsalida)
    salida.writerows(data)
    del salida
    csvsalida.close()

def getPolygon(distance, londiff, d, polar_positions, times, i):
    
    special = False
    arg = londiff * math.pi / (distance * 180)
    angle = math.asin(arg)  * 180 / math.pi

    pol_vertex = []
    for ind in range(2):
        p1 = [polar_positions[times[i]][1], polar_positions[times[i]][0]]
        vert1 = gcc.point_given_start_and_bearing(p1, angle - 90, d * 1000)
        vert2 = gcc.point_given_start_and_bearing(p1, angle + 90, d * 1000)            
        
        if vert1[0] < 0 and vert2[0] > 0:
            special = True
        elif vert2[0] < 0 and vert1[0] > 0:
            special = True

        vert1 = [vert1[1], vert1[0]]
        vert2 = [vert2[1], vert2[0]]

        pol_vertex.append(vert1)
        pol_vertex.append(vert2)

        i = i + 1
    
    r = LinearRing(pol_vertex)
    frame = Polygon(r)

    return frame, special, r

def sph2cartesian(latitude, longitude):
    R = 6371  # relative to centre of the earth
    longitude = longitude * math.pi / 180
    latitude = latitude * math.pi / 180

    X = R * math.cos(longitude) * math.cos(latitude)
    Y = R * math.sin(longitude) * math.cos(latitude)
    
    return X, Y

def Distance_objSats(time, time2, sat_position, objectives, objective):
    sat1 = np.array(sph2cartesian(sat_position[time][0], sat_position[time][1])) #lat, lon
    sat2 = np.array(sph2cartesian(sat_position[time2][0], sat_position[time2][1]))

    obj = np.array(sph2cartesian(objectives[objective][0], objectives[objective][1]))

    Line = [sat1, sat2]

    line = LineString(Line)
    p = Point(obj)

    return p.distance(line)

def saveObjective(objectivesInFootprint, objectives, objective, time, sat_position, percentage_distance):
    
    objectivesInFootprint[objective] = []

    for i in range(len(objectives[objective])):
        objectivesInFootprint[objective].append(objectives[objective][i])
    
    objectivesInFootprint[objective].append(percentage_distance)
    objectivesInFootprint[objective].append(sat_position[time])
    objectivesInFootprint[objective].append(time)
    
    return objectivesInFootprint

def objectivesInside(objectives, poly, objectivesInFootprint, time, time2, sat_position, special, corners, visualDistance):

    for objective in objectives:
        if objective not in objectivesInFootprint:
            node = Point(objectives[objective][0], objectives[objective][1])
            if poly.contains(node) and Distance_objSats(time, time2, sat_position, objectives, objective) <= visualDistance:
                percentage_distance = (visualDistance - Distance_objSats(time, time2, sat_position, objectives, objective)) * 100 / visualDistance

                if special:
                    limits = corners.bounds
                    if objectives[objective][0] > limits[0] and objectives[objective][0] < limits[2]:
                        if (limits[3] - limits[1]) < 180:
                            if objectives[objective][1] > limits[1] and objectives[objective][1] < limits[3]:
                                objectivesInFootprint = saveObjective(objectivesInFootprint, objectives, objective, time, sat_position, percentage_distance)
                        else:
                            if objectives[objective][1] < limits[1] and objectives[objective][1] > limits[3]:
                                objectivesInFootprint = saveObjective(objectivesInFootprint, objectives, objective, time, sat_position, percentage_distance)

                else:
                    objectivesInFootprint = saveObjective(objectivesInFootprint, objectives, objective, time, sat_position, percentage_distance)

    return objectivesInFootprint

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
    header.append("Time")

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
    salida = csv.writer(csvsalida)
    salida.writerows(data)
    del salida
    csvsalida.close()

def searchObjetives(times, polar_positions, visionDistance, objectivesInFootprint):

    print("Searching objetives under satellite footprint...")

    for i in range(len(times) - 1):
        print(times[i])
        lonDiff = polar_positions[times[i + 1]][1] - polar_positions[times[i]][1]

        #Arreglar cambio de cuadrantes
        if lonDiff < -180:
            lonDiff = abs(polar_positions[times[i + 1]][1] + polar_positions[times[i]][1])
        elif lonDiff > 180:
            lonDiff = - abs(polar_positions[times[i + 1]][1] + polar_positions[times[i]][1])

        lonDiff_distance = lonDiff * math.pi * Rtierra / 180

        p1 = [polar_positions[times[i]][1], polar_positions[times[i]][0]]
        p2 = [polar_positions[times[i + 1]][1], polar_positions[times[i + 1]][0]]
        distance = gcc.distance_between_points(p1, p2) / 1000

        #Obtener polígono de visión del sensor en el instante dado
        poly, special, corners = getPolygon(distance, lonDiff_distance, visionDistance, polar_positions, times, i)

        #Obtener ciudades en vista por el sensor en el instante dado
        objectivesInFootprint = objectivesInside(objectives, poly, objectivesInFootprint, times[i], times[i + 1], polar_positions, special, corners, visionDistance)

    return objectivesInFootprint

exportSat = init()

#GRAPHICAL USER INTERFACE
master = Tk()
master.title('Coverage calculation software')

actual_altitude = "Current altitude is " + str(alt) + ", introduce new altitude: "
altitude = Label(master, text = actual_altitude)
altitude.pack()

e1 = Entry(master)
e1.pack()

actual_inclination = "Current inclination is " + str(inc) + ", introduce new inclination: "
inclination = Label(master, text = actual_inclination)
inclination.pack()

e2 = Entry(master)
e2.pack()

actual_raan = "Current right ascension of the ascendint node is " + str(raan) + ", introduce new raan: "
raan2 = Label(master, text = actual_raan)
raan2.pack()

e3 = Entry(master)
e3.pack()

actual_FOV = "Current FOV is " + str(FOV) + " degrees, introduce new FOV: "
FOV2 = Label(master, text = actual_FOV)
FOV2.pack()

e4 = Entry(master)
e4.pack()

actual_cycle = "Current calculation period is " + str(cycle) + " hours, introduce new calculation period: "
cycle2 = Label(master, text = actual_cycle)
cycle2.pack()

e7 = Entry(master)
e7.pack()

if file == "":
    actual_fileinput = "Currently there's not file with objectives, introduce file name (example: worldcities.xlsx): "
else:
    actual_fileinput = "Current Excel with objectives is " + str(file) + ", introduce new file name: "

fileinput = Label(master, text = actual_fileinput)
fileinput.pack()

Label(master, text = "Make sure the file is in the same directory as the program").pack()

e5 = Entry(master)
e5.pack()

if objectives_init == "":
    actual_objectives = "Currently there're not objective's characteristics, introduce name, lat, lon and aditional characteristics: "
else:
    actual_objectives = "Current objective's characteristics are " + str(objectives_init) + ", introduce new objective's characteristics: "

objectives2 = Label(master, text = actual_objectives)
objectives2.pack()

Label(master, text = "Objective's characteristics should be the column's title from \ninput file and separated with commas, example: city_name,lat,lng,country").pack()

e6 = Entry(master)
e6.pack()

Label(master, text = "\n").pack()

button_changeparam = Button(master, text='Change parameters', width=25, command=get_parameters)
button_changeparam.pack()

button_exe = Button(master, text='EXECUTE', width=25, command=master.destroy, fg="green")
button_exe.pack()

button_exit = Button(master, text='EXIT', width=25, command=lambda: close(master), fg="red")
button_exit.pack()

Label(master, text = "\n").pack()

mainloop()

alt, inc, raan, FOV, file, objectives_init, cycle = checkValues(alt, inc, raan, FOV, file, objectives_init, cycle)

polar_positions = {}
objectivesInFootprint = {}

objectives = getObjectives(file, objectives_init)

visionDistance = visionArc(Rtierra, alt, FOV)
print(visionDistance)
exit()
#Creamos la órbita del satélite final con hora correcta
sat = Orbit.circular(Earth, alt, inc, raan)

#Propagamos el satélite en el tiempo
times = orbitTimes(cycle, timeDelta)

polar_positions = propagateSat(times, sat)

#Exportamos a un excel los datos de posición del satélite
if(exportSat):
    exportSatPositions(polar_positions)

objectivesInFootprint = searchObjetives(times, polar_positions, visionDistance, objectivesInFootprint) #igual falta meter objectives¿?

exportObjectives(objectivesInFootprint, objectives_init)
