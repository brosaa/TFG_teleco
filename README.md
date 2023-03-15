# TFG_teleco

Herramienta software de cálculo de cobertura de un sensor desde satélite en una órbita específica y localización de objetivos en la superficie terrestre

Autoría: Belén Rosa Alonso

Grado: Doble grado en ingeniería en tecnologías de las telecomunicaciones e ingeniería aeroespacial en aeronavegación

Software realizado con Thales Alenia Space como trabajo de fin de grado de ingeniería en tecnologías de las telecomunicaciones

Requisito de la instalación de los siguientes paquetes python:
- tkinter
- astropy
- numpy
- poliastro
- math
- pandas
- great_circle_calculator
- shapely
- csv

Código ejecutado con versión python 3.9

Antes de ejecutar, situar el terminal en los directorios calculodecobertura_capa_exportacionDeDatos\, calculodecobertura_capa_logica\ y calculodecobertura_capa_presentacion\ y ejecutar el comando  "pip install -e ." en cada uno de ellos.

Programa principal a ejecutar: calculo_coberturas_software.py, los parámetros de entrada se introducen en la interfaz gráfica posterior. 

## Docker
Para construir la imagen ejecutar en un terminal linux desde el home del proyecto:
   $docker build -t calculo-coberturas-software:latest -f Dockerfile .

Para ejecutar la imagen contruida anteriormente:
   $docker run -u=$(id -u $USER):$(id -g $USER) -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix:rw -it calculo-coberturas-software:latest



