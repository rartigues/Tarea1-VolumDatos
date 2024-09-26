# Tarea 1 - HyperLogLog y Jaccard

## Requisitos previos

- Compilador C++ que soporte C++11 o superior
- CMake (para compilar smhasher)

## Configuración

1. Clonar smhasher:
   ```
   git clone https://github.com/aappleby/smhasher.git
   ```

2. Copiar `MurmurHash3.cpp` y `MurmurHash3.h` al directorio de tu proyecto.

3. Colocar los archivos de genoma en el directorio del proyecto con los nombres:
   - genome1.txt
   - genome2.txt
   - genome3.txt
   - genome4.txt
   - genome5.txt

## Compilación

Compilar el programa principal:

```
g++ -std=c++11 -O3 main.cpp MurmurHash3.cpp -o hyperloglog_jaccard
```

## Ejecución

Ejecutar el programa:

```
./hyperloglog_jaccard
```

El programa procesará los genomas y mostrará los resultados de similitud y rendimiento para diferentes configuraciones de parámetros.

## Notas

- Asegúrate de que los archivos de genoma estén en formato de texto plano.
- El programa utiliza diferentes valores de k, w y b para la evaluación. Puedes modificar estos valores en el código fuente si es necesario.
