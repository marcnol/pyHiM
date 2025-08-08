#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classes and functions for file management
"""

import fcntl
import json
import multiprocessing
import os
import pathlib

import numpy as np
from dask.distributed import Client, LocalCluster, get_client

from core.pyhim_logging import print_log


class _GPUAssignPlugin:
    """Assigne un GPU unique par worker (round-robin) via un lock fichier."""

    def __init__(
        self, gpu_ids_env="PYHIM_GPU_IDS", state_path="/tmp/pyhim_gpu_rr.json"
    ):
        self.gpu_ids_env = gpu_ids_env
        self.state_path = state_path

    def setup(self, worker):
        gpu_ids = os.environ.get(self.gpu_ids_env, "0").split(",")
        gpu_ids = [g.strip() for g in gpu_ids if g.strip() != ""]
        if not gpu_ids:
            return  # CPU fallback

        state_file = pathlib.Path(self.state_path)
        state_file.parent.mkdir(parents=True, exist_ok=True)

        # section critique: réserver un ID GPU
        with open(state_file, "a+b") as f:
            f.seek(0)
            try:
                fcntl.flock(f, fcntl.LOCK_EX)
                try:
                    data = f.read().decode() or "{}"
                    state = json.loads(data)
                except Exception:
                    state = {}
                idx = state.get("idx", -1)
                next_idx = (idx + 1) % len(gpu_ids)
                chosen = gpu_ids[next_idx]
                state["idx"] = next_idx
                f.seek(0)
                f.truncate()
                f.write(json.dumps(state).encode())
                f.flush()
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)

        # **Très important** : fixer CUDA_VISIBLE_DEVICES avant tout import TF
        os.environ["CUDA_VISIBLE_DEVICES"] = str(chosen)
        os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"

        # log utile pour debug
        try:
            wid = getattr(worker, "id", "?")
            addr = getattr(worker, "address", "?")
        except Exception:
            wid, addr = "?", "?"
        print(f"[GPU-ASSIGN] worker={wid} addr={addr} -> CUDA_VISIBLE_DEVICES={chosen}")


class DaskCluster:
    """Used to manage parallel run thanks the Dask package"""

    def __init__(self, requested_nb_nodes, maximum_load=0.6, memory_per_worker=12000):
        self.requested_nb_nodes = requested_nb_nodes
        # self.n_threads will be created after exetution of initialize_cluster()
        self.n_threads = None
        self.maximum_load = maximum_load  # max number of workers that I can take
        self.memory_per_worker = memory_per_worker  # in Mb
        self.initialize_cluster()
        self.cluster = None
        self.client = None

    def initialize_cluster(self):
        """Defines the number of threads allocated"""
        number_cores_available = multiprocessing.cpu_count()

        # we want at least 12 GB per worker
        free_m = int(os.popen("free -t -m").readlines()[1].split()[-1])

        max_number_threads = int(
            np.min(
                [
                    number_cores_available * self.maximum_load,
                    free_m / self.memory_per_worker,
                ]
            )
        )

        self.n_threads = int(np.min([max_number_threads, self.requested_nb_nodes]))

        print_log(
            f"$ Cluster with {self.n_threads} workers started ({self.requested_nb_nodes} requested)"
        )

    def create_distributed_client(self):
        """Instance workers"""
        client = try_get_client()
        if client is not None:
            print_log("# Shutting down existing cluster! ")
            client.shutdown()
        else:
            print_log("$ No running cluster detected. Will start one.")

        self.cluster = LocalCluster(
            n_workers=self.n_threads,
            threads_per_worker=1,
            memory_limit="64GB",
            processes=True,
        )
        self.client = Client(self.cluster)
        # Enregistre le plugin sur chaque worker
        self.client.register_worker_plugin(_GPUAssignPlugin(), name="gpu-assign")
        print_log("$ Go to http://localhost:8787/status for information on progress...")


def try_get_client():
    """Check if client is alive

    Returns
    -------
    Dask.Client
        Client instance or None
    """
    try:
        client = get_client()
        client.restart()
    except ValueError:
        client = None

    return client
