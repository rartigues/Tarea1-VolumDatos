# Tarea 1 - Sketches para Estimación de Similitud entre Genomas

Este proyecto implementa y analiza algoritmos de streaming y sketches para la estimación de similitud entre genomas usando sketches de cardinalidad, específicamente HyperLogLog.

## Requisitos previos

- Compilador C++ que soporte C++11 o superior
- Git (para clonar el repositorio)

## Instalación

1. Clonar el repositorio:
   ```
   git clone https://github.com/rartigues/Tarea1-VolumDatos.git
   cd Tarea1-VolumDatos
   ```

2. Todos los archivos necesarios, incluyendo `MurmurHash3.cpp` y `MurmurHash3.h`, ya están incluidos en el repositorio.

## Compilación

Para compilar el proyecto, simplemente ejecuta el siguiente comando en el directorio del proyecto:

```
g++ -std=c++11 -O3 main.cpp MurmurHash3.cpp -o main
```

## Ejecución

Después de compilar, puedes ejecutar el programa con:

```
./main
```
También se puede ejecutar con el parametro `--log` para guardar los resultados en un archivo .csv

El programa procesará los genomas y mostrará los resultados de similitud y rendimiento para diferentes configuraciones de parámetros.

## Estructura del proyecto

- `main.cpp`: Contiene el código principal del proyecto.
- `MurmurHash3.cpp` y `MurmurHash3.h`: Implementación de la función hash MurmurHash3.
- `instances/`: Directorio que contiene los archivos de genoma para procesar.

## Notas

- El programa está configurado para procesar los genomas encontrados en el directorio `instances/`.
- Se utilizan diferentes valores de k, w y b para la evaluación, que pueden ser modificados en el código fuente si es necesario.
- Los resultados incluyen estimaciones de similitud de Jaccard usando k-mers y minimizers, así como métricas de error y tiempos de ejecución.

## Autores

- Roberto Artigues
- Oscar Castillo
