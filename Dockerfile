# Container for every step except 3D prediction (RhoFold+ has its own environment).
# Build: docker build -t aptamer-curation-pipeline:2.0 .
FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt
