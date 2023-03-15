FROM python:3.10-slim

RUN apt-get update -y && apt-get install tk -y

ADD dependencies  /tmp/dependencies
RUN python3 -m pip install --upgrade pip &&\
    pip3 install -r /tmp/dependencies/requirements.txt &&\
    pip3 install -e /tmp/dependencies/calculodecobertura_capa_exportacionDeDatos &&\
    pip3 install -e /tmp/dependencies/calculodecobertura_capa_logica &&\
    pip3 install -e /tmp/dependencies/calculodecobertura_capa_presentacion
    #rm -rf /tmp/dependencies

ADD app /app

CMD ["/app/calculo_coberturas_software.py"]
ENTRYPOINT ["python3"]

