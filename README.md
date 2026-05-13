# TFM — Reimplementación y Análisis del Algoritmo de Docking Proteína-Proteína ZDOCK en Python

**Autor:** belensalar  
**Universidad:** Universidad de Murcia  

---

## Descripción

Este repositorio contiene el código desarrollado en el Trabajo de Fin de Máster (TFM) sobre la reimplementación y análisis del algoritmo de docking rígido proteína-proteína ZDOCK en Python. El trabajo toma como referencia el artículo de Chen y Weng (2002) *"Docking Unbound Proteins Using Shape Complementarity, Desolvation, and Electrostatics"* y como punto de partida práctico una reimplementación del algoritmo en Julia desarrollada por Megan.

El objetivo principal es traducir el código de Julia a Python, aplicar estrategias de aceleración mediante el decorador `@njit` de Numba y el módulo `multiprocessing`, y evaluar el rendimiento del algoritmo sobre complejos proteína-proteína *unbound* de tipo anticuerpo-antígeno.

---

## Estructura del repositorio

```
📁 notebooks/
    📁 julia/
        proteins.jl          # Código original de Megan en Julia
    📁 python_original/
        proteins_1-3WD5.ipynb        # Ejecución 1 - Complejo 3WD5
        proteins_2-3WD5.ipynb        # Ejecución 2 - Complejo 3WD5
        proteins_3-3WD5.ipynb        # Ejecución 3 - Complejo 3WD5
        proteins_1-3MXW.ipynb        # Ejecución 1 - Complejo 3MXW
        proteins_2-3MXW.ipynb        # Ejecución 2 - Complejo 3MXW
        proteins_3-3MXW.ipynb        # Ejecución 3 - Complejo 3MXW
        proteins_1-5Y9J.ipynb        # Ejecución 1 - Complejo 5Y9J
        proteins_2-5Y9J.ipynb        # Ejecución 2 - Complejo 5Y9J
        proteins_3-5Y9J.ipynb        # Ejecución 3 - Complejo 5Y9J
        proteins_rand_3WD5.ipynb     # Ejecución aleatoria - Complejo 3WD5
        proteins_rand_3MXW.ipynb     # Ejecución aleatoria - Complejo 3MXW
        proteins_rand_5Y9J.ipynb     # Ejecución aleatoria - Complejo 5Y9J
    📁 optimizacion/
        cProfile_time.ipynb          # Análisis de rendimiento con cProfile
        Multiprocessing_Numba.ipynb  # Optimización con Numba y multiprocessing
📄 Dockerfile                        # Configuración del entorno Python reproducible
📄 README.md
```

---

## Entorno de ejecución

El código se ejecutó en el clúster HPC de la Universidad de Murcia, empleando el nodo de acceso **ibsen** y la cola de cómputo **Mendel**.

### Imagen Docker

La imagen Docker con todas las dependencias necesarias está disponible públicamente en Docker Hub:

```
docker pull belensalar/tfm-docking:v2
```

### Construcción del entorno desde el Dockerfile

Si se prefiere construir el entorno desde cero a partir del Dockerfile incluido en este repositorio:

```
docker build -t belensalar/tfm-docking:v2 .
```

### Ejecución del entorno

Para iniciar el contenedor y acceder al notebook de Jupyter:

```
docker run -it -p 8888:8888 --name tfm-notebook belensalar/tfm-docking:v2
```

Si se trabaja en un entorno HPC, es necesario configurar un túnel SSH para redirigir el puerto desde el nodo de cómputo al navegador local:

```
ssh -L 8888:<ip_nodo>:8888 <usuario>@<nodo_acceso>
```

Finalmente, abrir en el navegador:

```
http://localhost:8888
```

---

## Datos de entrada

Las estructuras proteicas empleadas corresponden a tres complejos anticuerpo-antígeno obtenidos del benchmark de anticuerpos de PierceLab:

- **3WD5**
- **3MXW**
- **5Y9J**

Los archivos PDB pueden descargarse desde el repositorio de PierceLab:
```
git clone https://github.com/piercelab/antibody_benchmark.git
```

---

## Dependencias principales

| Librería | Versión |
|----------|---------|
| Python | 3.11 |
| JupyterLab | 4.5.4 |
| NumPy | 2.4.2 |
| SciPy | 1.17.0 |
| pandas | 3.0.0 |
| Matplotlib | 3.10.8 |
| Numba | 0.64.0 |
| Biopython | 1.86 |

---

## Referencia

Chen, R., & Weng, Z. (2002). Docking unbound proteins using shape complementarity, desolvation, and electrostatics. *Proteins: Structure, Function, and Bioinformatics*, 47(3), 281-294.
