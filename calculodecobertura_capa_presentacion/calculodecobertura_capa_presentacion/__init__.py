"""calculodecobertura_capa_presentacion - Capa de presentación del software de cálculo de cobertura de un satélite, Trabajo de Fin de Grado en ingeniería de tecnologías de las telecomunicaciones URJC 2022 de Belén Rosa Alonso"""

__version__ = '0.0.1'
__author__ = 'Belén Rosa Alonso <b.rosaa.2017@alumnos.urjc.es>'
__all__ = []

def close(master):

    print("Closing program...")
    master.quit()
    exit()

def checkValues(alt, inc, raan, FOV, file, objectives_init, cycle):
    abort = False

    if(alt <= 0):
        print("Not valid altitude, please insert valid units for kilometers. \nClosing program...")
        abort = True    
    elif(file == ""):
        print("There is not an objetives file included. \nClosing program...")
        abort = True 
    elif(objectives_init == ""):
        print("There is not any objetive's characteristics to analyse included. \nClosing program...")
        abort = True 
    elif((cycle <= 0)):
        print("Not valid calculation time, please insert valid units for hours. \nClosing program...")
        abort = True
    
    if(abort):
        exit()
    
    return alt, inc, raan, FOV, file, objectives_init, cycle

def checkInput(i):

    if((i == "Y") or (i == "y")):
        exportSat = True
        return True, exportSat
    elif((i == "N") or (i == "n")):
        exportSat = False
        return True, exportSat
    else:
        return False, False

def init():

    i = input("Do you want to create an Excel file with the satellite positions [Y/N]: ")

    value, exportSat = checkInput(i)
    while(value == False):
        i = input("Incorrect input. \nDo you want to create an Excel file with the satellite positions [Y/N]: ")
        value, exportSat = checkInput(i)

    return exportSat

