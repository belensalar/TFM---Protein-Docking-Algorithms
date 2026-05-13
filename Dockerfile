# Imagen base
FROM python:3.11-slim

# Directorio de trabajo dentro del contenedor
WORKDIR /app

# Instalación de dependencias del sistema necesarias para compilación
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Instalación de librerías Python con versiones fijadas
RUN pip install --no-cache-dir \
    jupyterlab==4.5.4 \
    notebook==7.5.3 \
    numpy==2.4.2 \
    scipy==1.17.0 \
    pandas==3.0.0 \
    matplotlib==3.10.8 \
    numba==0.64.0 \
    llvmlite==0.46.0 \
    biopython==1.86 \
    plotly==6.5.2 \
    pyvista==0.47.0 \
    vtk==9.6.0

# Exponer el puerto de Jupyter
EXPOSE 8888

# Comando por defecto al iniciar el contenedor
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root"]
