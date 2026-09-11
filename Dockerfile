FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Variables d'environnement pour le support GPU / OpenGL (NVIDIA Container Toolkit)
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=graphics,utility,compute

# 1. Dépendances système, outils de build, bibliothèques OpenGL/GLX et X11
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    pkg-config \
    wget \
    sed \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    libgl1-mesa-glx \
    libgl1-mesa-dri \
    libglib2.0-0 \
    libsm6 \
    libxext6 \
    libxrender-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/CPlantBox

# 2. Copie du code source du dépôt
COPY . /opt/CPlantBox

# 3. Dépendances Python
RUN pip3 install --no-cache-dir --upgrade pip setuptools wheel && \
    pip3 install --no-cache-dir pybind11 && \
    pip3 install --no-cache-dir -r requirements.txt

# 4. Compilation et installation de CPlantBox
RUN cmake -B build -S . && \
    cmake --build build -j$(nproc) && \
    cmake --install build

# 5. Correctif automatique du bug d'API VTK (AddActor2D -> AddActor)
# S'exécute si le fichier est présent après l'installation
RUN if [ -f /usr/lib/python3/dist-packages/plantbox/visualisation/vtk_plot.py ]; then \
        sed -i 's/ren.AddActor2D(sb)/ren.AddActor(sb)/g' /usr/lib/python3/dist-packages/plantbox/visualisation/vtk_plot.py; \
    fi

ENV PYTHONPATH="/opt/CPlantBox:/opt/CPlantBox/src:${PYTHONPATH}"

CMD ["/bin/bash"]