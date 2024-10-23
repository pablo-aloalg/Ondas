[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gl/geoocean%2Fcourses%2Fwaves/HEAD?urlpath=tree/)

# ONDAS

Fernando J. Méndez Incera (fernando.mendez@unican.es)\
Jared Ortiz-Angulo Cantos (jared.ortizangulo@unican.es)

<a name="ins"></a>
## Data download
- - -
["DESCARGA DE DATOS"]

142 estaciones repartidas por el mundo (887MB) (serie temporal = 40 años; con datos horarios):

Link: https://e9kpv.app.goo.gl/48LnPegH1St7avYr5

Contraseña: geoocean

<a name="ins"></a>
## Install
- - -


### Activate conda environment


Descargamos todo el contenido de gitlab

Abrimos una terminal de conda: anconda prompt


Cambiamos al directorio en el que tenemos extraídos los archivos de gitlab

```
cd ondas-main
```

A partir del archivo environment.yml, creamos un nuevo env que se llama ondas2 (si da error, quitar el -f):

```
conda env create -f environment.yml
```


Activamos el environment:

```
conda activate ondas2

```

Lanzamos jupyter lab:

```
jupyter lab

```


