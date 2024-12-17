## Functions

import time
import sys
import h5py
import numpy as np
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau, CosineAnnealingLR, OneCycleLR

# Neural Operator Imports
from neuralop.models import FNO, TFNO
from neuralop import Trainer, LpLoss, H1Loss
from neuralop.datasets import load_darcy_flow_small
from neuralop.utils import count_model_params


def prepare_data_2D(h_all, p_train, p_val):
    h_all = np.real(h_all)

    # Input a
    a = h_all[:-1]

    # Target u (one time step after a)
    u = h_all[1:]
    n = len(h_all)

    # Split data into training, validation and test sets
    n_train = int(p_train * n)
    n_val = int(p_val * n)
    n_test = n - n_train - n_val

    print(f"n_train = {n_train}, n_val = {n_val}, n_test = {n_test}")

    train_x = torch.tensor(a[:n_train], dtype=torch.float32)  # Shape (n_train, 3, 256)
    val_x = torch.tensor(a[n_train:n_train + n_val], dtype=torch.float32)  # Shape (n_val, 3, 256)
    test_x = torch.tensor(a[n_train + n_val:], dtype=torch.float32)  # Shape (n_test, 3, 256)

    train_y = torch.tensor(u[:n_train], dtype=torch.float32)  # Shape (n_train, ...)
    val_y = torch.tensor(u[n_train:n_train + n_val], dtype=torch.float32)  # Shape (n_val, ...)
    test_y = torch.tensor(u[n_train + n_val:], dtype=torch.float32)  # Shape (n_test, ...)

    print(train_x.shape, train_y.shape, val_x.shape, val_y.shape, test_x.shape, test_y.shape)

    return n_train, n_val, n_test, a, u, n, train_x, train_y, val_x, val_y, test_x, test_y

def prepare_data_2D_FNO(h_all, p_train, p_val):
    h_all = np.real(h_all)

    # Input a
    a = h_all[:-1]
    a = a[:, np.newaxis, :, :]  # Add channel dimension

    # Target u (one time step after a)
    u = h_all[1:]
    u = u[:, np.newaxis, :, :]  # Add channel dimension

    n = len(h_all)
    
    # Split data into training, validation and test sets
    n_train = int(p_train * n)
    n_val = int(p_val * n)
    n_test = n - n_train - n_val

    print(f"n_train = {n_train}, n_val = {n_val}, n_test = {n_test}")

    train_x = torch.tensor(a[:n_train], dtype=torch.float32)  # Shape (n_train, 3, 256)
    val_x = torch.tensor(a[n_train:n_train + n_val], dtype=torch.float32)  # Shape (n_val, 3, 256)
    test_x = torch.tensor(a[n_train + n_val:], dtype=torch.float32)  # Shape (n_test, 3, 256)

    train_y = torch.tensor(u[:n_train], dtype=torch.float32)  # Shape (n_train, ...)
    val_y = torch.tensor(u[n_train:n_train + n_val], dtype=torch.float32)  # Shape (n_val, ...)
    test_y = torch.tensor(u[n_train + n_val:], dtype=torch.float32)  # Shape (n_test, ...)

    print(train_x.shape, train_y.shape, val_x.shape, val_y.shape, test_x.shape, test_y.shape)

    return n_train, n_val, n_test, a, u, n, train_x, train_y, val_x, val_y, test_x, test_y


# Prepare data in sequences
def prepare_data_sequences_2D(train_x, train_y, seq_length):
    # Fetch the dimensions of the input data
    timesteps, spatial_steps_x, spatial_steps_y = train_x.shape[0], train_x.shape[1], train_x.shape[2]

    # Number of sequences
    n_sequences = timesteps - seq_length

    # Prepare input and target sequences
    x_seq = torch.zeros((n_sequences, seq_length, spatial_steps_x, spatial_steps_y))
    y_seq = torch.zeros((n_sequences, spatial_steps_x, spatial_steps_y))

    for i in range(n_sequences):
        x_seq[i] = train_x[i:i + seq_length]
        y_seq[i] = train_y[i + seq_length]
    return x_seq, y_seq


def plot_losses(loss_train, loss_val, save, filename):
    plt.figure(figsize=(7, 4))
    plt.plot(loss_train, '.-', label='Training loss', linewidth=0.8, markersize=2)
    plt.plot(loss_val, '.-', label='Validation loss', linewidth=0.8, markersize=2)
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.title('Training and validation loss for each epoch')
    plt.xlabel('Epoch')
    plt.ylabel('Loss (MSE)')

    if save:
        plt.savefig(filename, format='pdf')

    plt.show()

def plot_error_2D(x, y, idx, t_all, pred_all, u, save, filename):
    # Prepare grid
    X, Y = np.meshgrid(x, y)
    cmap = 'Blues_r'

    fig = plt.figure(figsize=(12, 5))

    ax = fig.add_subplot(121, projection='3d')
    surf = ax.plot_surface(X, Y, np.abs(pred_all[idx]-u[idx]), cmap=cmap)
    ax.set_title(f'Absolute error (m) at t = {t_all[idx]:.2f}')
    ax.set_xlabel('Distance x (m)')
    ax.set_ylabel('Distance y (m)')
    ax.set_zlabel('Absolute error (m)')
    #fig.colorbar(surf, ax=ax, shrink=0.8, aspect=10, pad=0.15)  # Colorbar for 3D plot

    ax1 = fig.add_subplot(122)
    heatmap = ax1.imshow(np.abs(pred_all[idx]-u[idx]), cmap=cmap)
    ax1.set_title(f'Absolute error (m) at t = {t_all[idx]:.2f}')
    ax1.set_xlabel('Distance x (m)')
    ax1.set_ylabel('Distance y (m)')
    fig.colorbar(heatmap, ax=ax1, shrink=0.8, aspect=10, pad = 0.15) # Add colorbar

    plt.tight_layout()
    if save:
        plt.savefig(filename)

    plt.show()





