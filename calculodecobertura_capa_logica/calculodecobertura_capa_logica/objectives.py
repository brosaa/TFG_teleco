"""calculodecobertura_capa_logica - Capa lógica o de negocios del software de cálculo de cobertura de un satélite, Trabajo de Fin de Grado en ingeniería de tecnologías de las telecomunicaciones URJC 2022 de Belén Rosa Alonso"""

__version__ = '0.0.1'
__author__ = 'Belén Rosa Alonso <b.rosaa.2017@alumnos.urjc.es>'
__all__ = []

import pandas as pd
from math import pi, asin, cos, sin
import great_circle_calculator.great_circle_calculator as gcc
from shapely.geometry.polygon import Polygon, LinearRing, Point, LineString
#Probar a hacer visionarc desde aqui
import numpy as np

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

def getPolygon(distance, londiff, d, polar_positions, times, i):
    
    special = False
    arg = londiff * pi / (distance * 180)
    angle = asin(arg)  * 180 / pi

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
    longitude = longitude * pi / 180
    latitude = latitude * pi / 180

    X = R * cos(longitude) * cos(latitude)
    Y = R * sin(longitude) * cos(latitude)
    
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

def searchObjetives(times, polar_positions, visionDistance, objectivesInFootprint, Rtierra, file, objectives_init):

    print("Searching objetives under satellite footprint...")

    for i in range(len(times) - 1):
        if(times[i][len(times[i]) - 1] == "0"):
            print("Calculating at instant " + times[i])

        lonDiff = polar_positions[times[i + 1]][1] - polar_positions[times[i]][1]

        #Arreglar cambio de cuadrantes
        if lonDiff < -180:
            lonDiff = abs(polar_positions[times[i + 1]][1] + polar_positions[times[i]][1])
        elif lonDiff > 180:
            lonDiff = - abs(polar_positions[times[i + 1]][1] + polar_positions[times[i]][1])

        lonDiff_distance = lonDiff * pi * Rtierra / 180

        p1 = [polar_positions[times[i]][1], polar_positions[times[i]][0]]
        p2 = [polar_positions[times[i + 1]][1], polar_positions[times[i + 1]][0]]
        distance = gcc.distance_between_points(p1, p2) / 1000

        #Obtener polígono de visión del sensor en el instante dado
        poly, special, corners = getPolygon(distance, lonDiff_distance, visionDistance, polar_positions, times, i)

        objectives = getObjectives(file, objectives_init)

        #Obtener ciudades en vista por el sensor en el instante dado
        objectivesInFootprint = objectivesInside(objectives, poly, objectivesInFootprint, times[i], times[i + 1], polar_positions, special, corners, visionDistance)

    return objectivesInFootprint