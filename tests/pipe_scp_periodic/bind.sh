#!/bin/bash

# Make sure OMPI_COMM_WORLD_LOCAL_RANK is defined
if [ -z "$OMPI_COMM_WORLD_LOCAL_RANK" ]; then
  echo "[ERROR] OMPI_COMM_WORLD_LOCAL_RANK not set. Exiting."
  exit 1
fi

# Map MPI local rank to GPU ID 0–3
export LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
export CUDA_VISIBLE_DEVICES=${LOCAL_RANK}

echo "[LOG] local rank $LOCAL_RANK: bind to GPU $CUDA_VISIBLE_DEVICES"
echo ""

# Launch the actual command
exec "$@"

