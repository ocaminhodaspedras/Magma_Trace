"""
Magma_Trace: Simulador de Evolução Geoquímica de Elementos-Traço
================================================================
Autor: Felipe Padilha Leitzke
Afiliação: Departamento de Mineralogia e Petrologia (DEMIPE),
           Instituto de Geociências (IGEO), Universidade Federal do Rio Grande do Sul (UFRGS)
Contato para sugestões, dúvidas, erros, etc: felipe.leitzke@ufrgs.br
Ano: 2026
Licença: MIT License

Escopo Científico e Metodológico:
---------------------------------
Aplicativo desktop didático interativo para modelagem quantitativa do fracionamento de elementos-traço
(em razões relativas Cl/C0, concentrações absolutas em ppm e padrões normalizados ao Condrito C1)
durante processos ígneos de fusão parcial (modal, não-modal cotética e não-modal peritética,
em equilíbrio e fracionada; Shaw, 1970), cristalização fracionada de Rayleigh e assimilação
acoplada à cristalização fracionada (AFC; DePaolo, 1981).

Nota sobre o Banco de Dados (Kd, C0 e Ca):
------------------------------------------
A base interna de coeficientes de partição (DEFAULT_KD_DATABASE) constitui uma compilação
didática derivada de Arth (1976), Fujimaki et al. (1984), Rollinson (1993) e Rollinson & Pease (2021),
com interpolação na parábola de Onuma para aqueles sem determinação
experimental. Os reservatórios de concentração inicial C0 (ppm), assimilante Ca (ppm) e Condrito C1
seguem Sun & McDonough (1989, Tabela 1), McDonough & Sun (1995) e Rudnick & Gao (2003, Tabela 3).
Para aplicações estritas de pesquisa, recomenda-se a importação de dados específicos via arquivo CSV.

Declaração de Transparência no Uso de IA:
-----------------------------------------
A concepção petrológica, a definição das restrições termodinâmicas e de balanço de massa
(incluindo limites de exaustão de fases F_max e reações incongruentes), o tratamento de
singularidades numéricas e a auditoria matemática foram conduzidos integralmente pelo autor.
A implementação sintática da interface gráfica (PyQt5) e a estruturação do código em Python
contaram com assistência de modelos de linguagem de grande escala (LLMs), sob supervisão e
validação humana.
"""

import os
import sys
import copy
import struct
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QSlider, QGroupBox, QFormLayout,
                             QComboBox, QPushButton, QMessageBox, QDialog, QListWidget,
                             QListWidgetItem, QTableWidget, QTableWidgetItem, QHeaderView,
                             QTextBrowser, QFileDialog, QGridLayout, QCheckBox)
from PyQt5.QtCore import Qt, QByteArray, QBuffer, QIODevice, QRect
from PyQt5.QtGui import QPixmap, QPainter, QColor, QFont, QPen, QLinearGradient, QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.colors as mcolors


ELEMENTOS_MASTER = [
    'Rb', 'Ba', 'Sr',
    'Nb', 'Zr', 'Ti', 'Y',
    'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Sc', 'V', 'Cr', 'Ni'
]

LISTA_ETR = ['La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']

ORDEM_SPIDERGRAM = [
    'Rb', 'Ba', 'Nb', 'La', 'Ce', 'Pr', 'Sr', 'Nd', 'Zr', 'Sm', 'Eu', 'Ti',
    'Gd', 'Tb', 'Dy', 'Y', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Sc', 'V', 'Cr', 'Ni'
]

TODOS_MINERAIS = ['Olivina', 'Opx', 'Cpx', 'Granada', 'Plagioclásio', 'Anfibólio', 'Biotita']
CORES_HEX = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())[10:40]

MODA_PADRAO_SISTEMA = {
    'Basáltico (Máfico)':         {'Olivina': 55, 'Opx': 25, 'Cpx': 15, 'Granada': 5},
    'Andesítico (Intermediário)': {'Plagioclásio': 50, 'Cpx': 25, 'Opx': 15, 'Anfibólio': 10},
    'Riolítico (Félsico)':        {'Plagioclásio': 50, 'Anfibólio': 25, 'Biotita': 15, 'Cpx': 10}
}

# Condrito C1 em ppm: ETR, LILE e HFSE de Sun & McDonough (1989, Tabela 1);
# Metais de transição (Sc, V, Cr, Ni) de McDonough & Sun (1995, Tabela 1).
CONDRITO_C1_PPM = {
    'Rb': 2.32, 'Ba': 2.41, 'Sr': 7.26, 'Nb': 0.246, 'Zr': 3.87, 'Ti': 440.0, 'Y': 1.56,
    'La': 0.237, 'Ce': 0.612, 'Pr': 0.095, 'Nd': 0.467, 'Sm': 0.153, 'Eu': 0.058, 'Gd': 0.2055,
    'Tb': 0.0374, 'Dy': 0.254, 'Ho': 0.0566, 'Er': 0.1655, 'Tm': 0.0255, 'Yb': 0.170, 'Lu': 0.0254,
    'Sc': 5.92, 'V': 56.0, 'Cr': 2650.0, 'Ni': 10800.0
}

# Reservatórios geoquímicos de referência em ppm:
# Manto Primitivo, N-MORB e OIB de Sun & McDonough (1989, Tabela 1) com Sc, V, Cr, Ni de McDonough & Sun (1995)
# e médias de MORB/OIB; Crosta Continental Superior de Rudnick & Gao (2003, Tabela 3).
RESERVATORIOS_C0_PPM = {
    'Manto Primitivo (Sun & McDonough, 1989)': {
        'Rb': 0.635, 'Ba': 6.989, 'Sr': 21.1, 'Nb': 0.713, 'Zr': 11.2, 'Ti': 1300.0, 'Y': 4.55,
        'La': 0.687, 'Ce': 1.775, 'Pr': 0.276, 'Nd': 1.354, 'Sm': 0.444, 'Eu': 0.168, 'Gd': 0.596,
        'Tb': 0.108, 'Dy': 0.737, 'Ho': 0.164, 'Er': 0.480, 'Tm': 0.074, 'Yb': 0.493, 'Lu': 0.074,
        'Sc': 16.2, 'V': 82.0, 'Cr': 2625.0, 'Ni': 1960.0
    },
    'Basalto N-MORB (Sun & McDonough, 1989)': {
        'Rb': 0.56, 'Ba': 6.30, 'Sr': 90.0, 'Nb': 2.33, 'Zr': 74.0, 'Ti': 7600.0, 'Y': 28.0,
        'La': 2.50, 'Ce': 7.50, 'Pr': 1.32, 'Nd': 7.30, 'Sm': 2.63, 'Eu': 1.02, 'Gd': 3.68,
        'Tb': 0.67, 'Dy': 4.55, 'Ho': 1.01, 'Er': 2.97, 'Tm': 0.456, 'Yb': 3.05, 'Lu': 0.455,
        'Sc': 40.0, 'V': 250.0, 'Cr': 250.0, 'Ni': 100.0
    },
    'Basalto OIB (Sun & McDonough, 1989)': {
        'Rb': 31.0, 'Ba': 350.0, 'Sr': 660.0, 'Nb': 48.0, 'Zr': 280.0, 'Ti': 17200.0, 'Y': 29.0,
        'La': 37.0, 'Ce': 80.0, 'Pr': 9.70, 'Nd': 38.5, 'Sm': 10.0, 'Eu': 3.00, 'Gd': 7.62,
        'Tb': 1.05, 'Dy': 5.60, 'Ho': 1.06, 'Er': 2.62, 'Tm': 0.35, 'Yb': 2.16, 'Lu': 0.30,
        'Sc': 25.0, 'V': 220.0, 'Cr': 150.0, 'Ni': 120.0
    },
    'Crosta Continental Superior (Rudnick & Gao, 2003)': {
        'Rb': 84.0, 'Ba': 628.0, 'Sr': 320.0, 'Nb': 12.0, 'Zr': 193.0, 'Ti': 3800.0, 'Y': 21.0,
        'La': 31.0, 'Ce': 63.0, 'Pr': 7.10, 'Nd': 27.0, 'Sm': 4.70, 'Eu': 1.00, 'Gd': 4.00,
        'Tb': 0.70, 'Dy': 3.90, 'Ho': 0.83, 'Er': 2.30, 'Tm': 0.30, 'Yb': 1.96, 'Lu': 0.31,
        'Sc': 14.0, 'V': 97.0, 'Cr': 92.0, 'Ni': 47.0
    }
}

DEFAULT_KD_DATABASE = {
    'Basáltico (Máfico)': {
        'Olivina':      {'Rb': 0.00018, 'Ba': 0.00010, 'Sr': 0.00019, 'Nb': 0.0013, 'Zr': 0.0080, 'Ti': 0.0096, 'Y': 0.0068, 'La': 0.000022, 'Ce': 0.000045, 'Pr': 0.000085, 'Nd': 0.00020, 'Sm': 0.00064, 'Eu': 0.00075, 'Gd': 0.0012, 'Tb': 0.00220, 'Dy': 0.0040, 'Ho': 0.00643, 'Er': 0.0110, 'Tm': 0.0140, 'Yb': 0.0188, 'Lu': 0.0280, 'Sc': 0.12, 'V': 0.0896, 'Cr': 1.18, 'Ni': 24.01},
        'Opx':          {'Rb': 0.0030, 'Ba': 0.000006, 'Sr': 0.00120, 'Nb': 0.0040, 'Zr': 0.0320, 'Ti': 0.1400, 'Y': 0.0950, 'La': 0.000800, 'Ce': 0.001600, 'Pr': 0.003200, 'Nd': 0.00560, 'Sm': 0.01500, 'Eu': 0.02200, 'Gd': 0.0340, 'Tb': 0.05400, 'Dy': 0.0770, 'Ho': 0.10000, 'Er': 0.1200, 'Tm': 0.1700, 'Yb': 0.2200, 'Lu': 0.2200, 'Sc': 1.29, 'V': 0.8560, 'Cr': 3.52, 'Ni': 7.38},
        'Cpx':          {'Rb': 0.0100, 'Ba': 0.000200, 'Sr': 0.08800, 'Nb': 0.0037, 'Zr': 0.1110, 'Ti': 0.3800, 'Y': 1.3500, 'La': 0.086000, 'Ce': 0.175000, 'Pr': 0.289000, 'Nd': 0.47000, 'Sm': 0.81000, 'Eu': 0.91000, 'Gd': 1.0400, 'Tb': 1.22000, 'Dy': 1.4000, 'Ho': 1.35000, 'Er': 1.3400, 'Tm': 1.3800, 'Yb': 1.4200, 'Lu': 1.1600, 'Sc': 1.75, 'V': 2.9000, 'Cr': 8.10, 'Ni': 2.60},
        'Granada':      {'Rb': 0.0007, 'Ba': 0.000370, 'Sr': 0.00740, 'Nb': 0.0055, 'Zr': 0.7300, 'Ti': 0.2900, 'Y': 8.3000, 'La': 0.004700, 'Ce': 0.017900, 'Pr': 0.059300, 'Nd': 0.17000, 'Sm': 0.87000, 'Eu': 1.63000, 'Gd': 2.5500, 'Tb': 4.20000, 'Dy': 6.2000, 'Ho': 8.20000, 'Er': 9.6000, 'Tm': 11.100, 'Yb': 12.600, 'Lu': 13.700, 'Sc': 2.62, 'V': 3.6000, 'Cr': 3.80, 'Ni': 0.60},
        'Plagioclásio': {'Rb': 0.0710, 'Ba': 0.230000, 'Sr': 1.83000, 'Nb': 0.0100, 'Zr': 0.0010, 'Ti': 0.0380, 'Y': 0.0260, 'La': 0.148000, 'Ce': 0.111000, 'Pr': 0.095000, 'Nd': 0.08100, 'Sm': 0.06700, 'Eu': 0.34000, 'Gd': 0.0630, 'Tb': 0.04500, 'Dy': 0.0290, 'Ho': 0.02500, 'Er': 0.0180, 'Tm': 0.0150, 'Yb': 0.0110, 'Lu': 0.0090, 'Sc': 0.04, 'V': 0.0100, 'Cr': 0.04, 'Ni': 0.04},
        'Anfibólio':    {'Rb': 0.0900, 'Ba': 0.385000, 'Sr': 0.46000, 'Nb': 0.3900, 'Zr': 0.3700, 'Ti': 2.6650, 'Y': 1.3600, 'La': 0.200000, 'Ce': 0.350000, 'Pr': 0.510000, 'Nd': 0.73000, 'Sm': 1.07500, 'Eu': 1.13000, 'Gd': 1.3700, 'Tb': 1.39000, 'Dy': 1.4050, 'Ho': 1.36000, 'Er': 1.2750, 'Tm': 1.2000, 'Yb': 1.0500, 'Lu': 0.8500, 'Sc': 4.94, 'V': 6.0800, 'Cr': 6.03, 'Ni': 6.80}
    },
    'Andesítico (Intermediário)': {
        'Olivina':      {'Rb': 0.010, 'Ba': 0.010, 'Sr': 0.010, 'Nb': 0.010, 'Zr': 0.010, 'Ti': 0.030, 'Y': 0.038, 'La': 0.007, 'Ce': 0.008, 'Pr': 0.010, 'Nd': 0.011, 'Sm': 0.015, 'Eu': 0.016, 'Gd': 0.019, 'Tb': 0.025, 'Dy': 0.030, 'Ho': 0.037, 'Er': 0.044, 'Tm': 0.056, 'Yb': 0.070, 'Lu': 0.110, 'Sc': 0.28, 'V': 0.080, 'Cr': 2.80, 'Ni': 20.8},
        'Opx':          {'Rb': 0.022, 'Ba': 0.013, 'Sr': 0.010, 'Nb': 0.027, 'Zr': 0.031, 'Ti': 0.250, 'Y': 0.340, 'La': 0.002, 'Ce': 0.005, 'Pr': 0.010, 'Nd': 0.019, 'Sm': 0.050, 'Eu': 0.057, 'Gd': 0.140, 'Tb': 0.200, 'Dy': 0.260, 'Ho': 0.340, 'Er': 0.430, 'Tm': 0.510, 'Yb': 0.590, 'Lu': 0.750, 'Sc': 3.87, 'V': 1.000, 'Cr': 13.0, 'Ni': 12.0},
        'Cpx':          {'Rb': 0.030, 'Ba': 0.020, 'Sr': 0.280, 'Nb': 0.025, 'Zr': 0.290, 'Ti': 0.400, 'Y': 1.350, 'La': 0.099, 'Ce': 0.193, 'Pr': 0.290, 'Nd': 0.420, 'Sm': 0.750, 'Eu': 0.820, 'Gd': 0.910, 'Tb': 1.050, 'Dy': 1.200, 'Ho': 1.350, 'Er': 1.300, 'Tm': 1.150, 'Yb': 0.900, 'Lu': 0.800, 'Sc': 3.20, 'V': 3.150, 'Cr': 34.0, 'Ni': 6.80},
        'Granada':      {'Rb': 0.010, 'Ba': 0.010, 'Sr': 0.015, 'Nb': 0.040, 'Zr': 0.530, 'Ti': 1.200, 'Y': 14.50, 'La': 0.028, 'Ce': 0.080, 'Pr': 0.130, 'Nd': 0.222, 'Sm': 0.810, 'Eu': 1.540, 'Gd': 4.590, 'Tb': 6.300, 'Dy': 9.000, 'Ho': 15.00, 'Er': 20.00, 'Tm': 22.00, 'Yb': 24.00, 'Lu': 24.00, 'Sc': 7.00, 'V': 8.000, 'Cr': 22.0, 'Ni': 0.60},
        'Plagioclásio': {'Rb': 0.070, 'Ba': 0.360, 'Sr': 2.840, 'Nb': 0.025, 'Zr': 0.010, 'Ti': 0.047, 'Y': 0.019, 'La': 0.228, 'Ce': 0.170, 'Pr': 0.130, 'Nd': 0.100, 'Sm': 0.065, 'Eu': 1.210, 'Gd': 0.052, 'Tb': 0.038, 'Dy': 0.025, 'Ho': 0.018, 'Er': 0.014, 'Tm': 0.012, 'Yb': 0.010, 'Lu': 0.009, 'Sc': 0.01, 'V': 0.050, 'Cr': 0.02, 'Ni': 0.01},
        'Anfibólio':    {'Rb': 0.140, 'Ba': 0.120, 'Sr': 0.280, 'Nb': 0.800, 'Zr': 0.400, 'Ti': 3.000, 'Y': 1.650, 'La': 0.170, 'Ce': 0.260, 'Pr': 0.420, 'Nd': 0.630, 'Sm': 1.370, 'Eu': 1.210, 'Gd': 1.490, 'Tb': 1.650, 'Dy': 1.770, 'Ho': 1.600, 'Er': 1.470, 'Tm': 1.300, 'Yb': 1.150, 'Lu': 1.070, 'Sc': 10.5, 'V': 9.650, 'Cr': 30.0, 'Ni': 10.0}
    },
    'Riolítico (Félsico)': {
        'Opx':          {'Rb': 0.003, 'Ba': 0.003, 'Sr': 0.009, 'Nb': 0.150, 'Zr': 0.120, 'Ti': 0.400, 'Y': 0.480, 'La': 0.030, 'Ce': 0.050, 'Pr': 0.070, 'Nd': 0.100, 'Sm': 0.180, 'Eu': 0.150, 'Gd': 0.250, 'Tb': 0.320, 'Dy': 0.400, 'Ho': 0.500, 'Er': 0.620, 'Tm': 0.720, 'Yb': 0.850, 'Lu': 0.900, 'Sc': 8.50, 'V': 2.600, 'Cr': 10.0, 'Ni': 8.00},
        'Cpx':          {'Rb': 0.032, 'Ba': 0.131, 'Sr': 0.516, 'Nb': 0.300, 'Zr': 0.500, 'Ti': 0.700, 'Y': 3.650, 'La': 0.515, 'Ce': 0.646, 'Pr': 0.900, 'Nd': 1.400, 'Sm': 2.700, 'Eu': 1.560, 'Gd': 3.100, 'Tb': 3.400, 'Dy': 3.700, 'Ho': 3.600, 'Er': 3.500, 'Tm': 3.200, 'Yb': 2.800, 'Lu': 2.700, 'Sc': 22.0, 'V': 7.440, 'Cr': 30.0, 'Ni': 4.50},
        'Granada':      {'Rb': 0.009, 'Ba': 0.015, 'Sr': 0.015, 'Nb': 0.100, 'Zr': 1.200, 'Ti': 1.200, 'Y': 33.00, 'La': 0.070, 'Ce': 0.150, 'Pr': 0.300, 'Nd': 0.520, 'Sm': 2.660, 'Eu': 1.500, 'Gd': 10.50, 'Tb': 18.00, 'Dy': 28.60, 'Ho': 33.00, 'Er': 38.00, 'Tm': 40.00, 'Yb': 43.00, 'Lu': 40.00, 'Sc': 16.0, 'V': 7.000, 'Cr': 12.0, 'Ni': 0.60},
        'Plagioclásio': {'Rb': 0.041, 'Ba': 0.308, 'Sr': 4.400, 'Nb': 0.060, 'Zr': 0.048, 'Ti': 0.050, 'Y': 0.060, 'La': 0.380, 'Ce': 0.270, 'Pr': 0.230, 'Nd': 0.190, 'Sm': 0.130, 'Eu': 2.150, 'Gd': 0.097, 'Tb': 0.075, 'Dy': 0.064, 'Ho': 0.058, 'Er': 0.055, 'Tm': 0.050, 'Yb': 0.049, 'Lu': 0.046, 'Sc': 0.04, 'V': 0.010, 'Cr': 0.02, 'Ni': 0.01},
        'Anfibólio':    {'Rb': 0.045, 'Ba': 0.100, 'Sr': 0.120, 'Nb': 1.300, 'Zr': 1.200, 'Ti': 7.000, 'Y': 7.600, 'La': 0.540, 'Ce': 0.900, 'Pr': 1.600, 'Nd': 2.800, 'Sm': 5.400, 'Eu': 3.600, 'Gd': 7.000, 'Tb': 7.500, 'Dy': 8.200, 'Ho': 7.800, 'Er': 7.000, 'Tm': 6.000, 'Yb': 5.000, 'Lu': 4.500, 'Sc': 20.0, 'V': 12.00, 'Cr': 30.0, 'Ni': 8.00},
        'Biotita':      {'Rb': 3.260, 'Ba': 6.360, 'Sr': 0.120, 'Nb': 4.600, 'Zr': 1.200, 'Ti': 15.00, 'Y': 0.320, 'La': 0.320, 'Ce': 0.320, 'Pr': 0.310, 'Nd': 0.290, 'Sm': 0.260, 'Eu': 0.240, 'Gd': 0.280, 'Tb': 0.290, 'Dy': 0.290, 'Ho': 0.320, 'Er': 0.350, 'Tm': 0.400, 'Yb': 0.440, 'Lu': 0.330, 'Sc': 11.3, 'V': 18.00, 'Cr': 12.6, 'Ni': 3.50}
    }
}

KD_DATABASE = copy.deepcopy(DEFAULT_KD_DATABASE)
C0_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM['Manto Primitivo (Sun & McDonough, 1989)'])
CA_PPM_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM['Crosta Continental Superior (Rudnick & Gao, 2003)'])


def construir_icone_mt(caminho_ico="magma_trace.ico"):
    """
    Desenha o emblema 'M-T' em resolução 256x256, retorna um objeto QIcon
    para a janela em execução e grava um arquivo binário .ico no disco
    (apenas quando executado como script .py) para compilação com PyInstaller.
    """
    tamanho = 256
    pixmap = QPixmap(tamanho, tamanho)
    pixmap.fill(Qt.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.Antialiasing, True)
    painter.setRenderHint(QPainter.TextAntialiasing, True)

    grad = QLinearGradient(0, 0, tamanho, tamanho)
    grad.setColorAt(0.0, QColor("#1c1c1c"))
    grad.setColorAt(1.0, QColor("#2d1410"))

    pen_borda = QPen(QColor("#b8331a"))
    pen_borda.setWidth(14)
    painter.setPen(pen_borda)
    painter.setBrush(grad)
    painter.drawRoundedRect(10, 10, tamanho - 20, tamanho - 20, 36, 36)

    pen_linha = QPen(QColor("#ff5733"))
    pen_linha.setWidth(5)
    painter.setPen(pen_linha)
    painter.drawLine(36, 198, tamanho - 36, 198)

    fonte = QFont("Segoe UI", 68, QFont.Bold)
    painter.setFont(fonte)
    painter.setPen(QColor("#f2f2f2"))
    painter.drawText(QRect(0, 22, tamanho, 165), Qt.AlignCenter, "M-T")

    painter.end()

    if not getattr(sys, 'frozen', False) and not os.path.exists(caminho_ico):
        try:
            ba = QByteArray()
            buf = QBuffer(ba)
            buf.open(QIODevice.WriteOnly)
            pixmap.save(buf, "PNG")
            png_bytes = bytes(ba)

            icondir = struct.pack('<HHH', 0, 1, 1)
            icondirentry = struct.pack('<BBBBHHII', 0, 0, 0, 0, 1, 32, len(png_bytes), 22)
            with open(caminho_ico, 'wb') as f_ico:
                f_ico.write(icondir + icondirentry + png_bytes)
        except Exception:
            pass

    return QIcon(pixmap), pixmap


class TabelaKdDialog(QDialog):
    def __init__(self, callback_atualizar=None, callback_c0_custom=None, callback_ca_custom=None, callback_reset=None):
        super().__init__()
        self.callback_atualizar = callback_atualizar
        self.callback_c0_custom = callback_c0_custom
        self.callback_ca_custom = callback_ca_custom
        self.callback_reset = callback_reset
        self.setWindowTitle("Parâmetros Geoquímicos: Kd, Fonte/Magma Inicial C0 (ppm) e Assimilante Ca (ppm)")
        self.resize(1160, 610)
        layout = QVBoxLayout()

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel("<b>Sistema Magmático:</b>"))
        self.cb_magma_tab = QComboBox()
        self.cb_magma_tab.addItems(list(KD_DATABASE.keys()))
        self.cb_magma_tab.currentIndexChanged.connect(self.carregar_tabela)
        top_layout.addWidget(self.cb_magma_tab, 1)
        layout.addLayout(top_layout)

        self.tabela = QTableWidget()
        self.tabela.itemChanged.connect(self.celula_editada)
        layout.addWidget(self.tabela)

        lbl_nota = QLabel(
            "<small><i>Nota Metodológica: Os coeficientes de partição (Kd) constituem uma compilação didática derivada de Arth (1976), "
            "Fujimaki et al. (1984), Rollinson (1993) e Rollinson & Pease (2021). 'C0 (ppm)' define a concentração inicial na fonte ou magma parental "
            "(padrão: Manto Primitivo de Sun & McDonough, 1989) e 'Ca (ppm)' define a concentração na rocha encaixante assimilada no modelo AFC "
            "(padrão: Crosta Continental Superior de Rudnick & Gao, 2003). A coluna 'Ca/C0 (Calc.)' exibe a razão resultante automaticamente.</i></small>"
        )
        lbl_nota.setWordWrap(True)
        layout.addWidget(lbl_nota)

        botoes_layout = QHBoxLayout()
        btn_importar = QPushButton("Importar CSV (Padrão / Excel PT-BR)")
        btn_importar.clicked.connect(self.importar_csv)
        btn_exportar = QPushButton("Exportar CSV (Sistema Atual)")
        btn_exportar.clicked.connect(self.exportar_csv)
        btn_reset = QPushButton("Restaurar Padrão")
        btn_reset.clicked.connect(self.restaurar_padrao)
        btn_fechar = QPushButton("Aplicar e Fechar")
        btn_fechar.clicked.connect(self.accept)

        botoes_layout.addWidget(btn_importar)
        botoes_layout.addWidget(btn_exportar)
        botoes_layout.addWidget(btn_reset)
        botoes_layout.addStretch()
        botoes_layout.addWidget(btn_fechar)
        layout.addLayout(botoes_layout)

        self.setLayout(layout)
        self.carregar_tabela()

    def carregar_tabela(self):
        self.tabela.blockSignals(True)
        magma = self.cb_magma_tab.currentText()
        minerais = list(KD_DATABASE[magma].keys())
        colunas = minerais + ["C0 (ppm)", "Ca (ppm)", "Ca/C0 (Calc.)"]
        self.tabela.setRowCount(len(ELEMENTOS_MASTER))
        self.tabela.setColumnCount(len(colunas))
        self.tabela.setHorizontalHeaderLabels(colunas)
        self.tabela.setVerticalHeaderLabels(ELEMENTOS_MASTER)

        for row, el in enumerate(ELEMENTOS_MASTER):
            for col, mineral in enumerate(minerais):
                val = KD_DATABASE[magma][mineral].get(el, 0.0)
                item = QTableWidgetItem(f"{val:.5f}")
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.tabela.setItem(row, col, item)

            val_c0 = C0_DATABASE.get(el, 1.0)
            item_c0 = QTableWidgetItem(f"{val_c0:.4g}")
            item_c0.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.tabela.setItem(row, len(minerais), item_c0)

            val_ca = CA_PPM_DATABASE.get(el, 1.0)
            item_ca = QTableWidgetItem(f"{val_ca:.4g}")
            item_ca.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.tabela.setItem(row, len(minerais) + 1, item_ca)

            razao_ca_c0 = val_ca / max(1e-12, val_c0)
            item_raz = QTableWidgetItem(f"{razao_ca_c0:.3f}")
            item_raz.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            item_raz.setFlags(item_raz.flags() & ~Qt.ItemIsEditable)
            item_raz.setBackground(QColor("#f0f0f0"))
            self.tabela.setItem(row, len(minerais) + 2, item_raz)

        self.tabela.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.tabela.blockSignals(False)

    def celula_editada(self, item):
        magma = self.cb_magma_tab.currentText()
        minerais = list(KD_DATABASE[magma].keys())
        row = item.row()
        col = item.column()
        el = ELEMENTOS_MASTER[row]

        if col == len(minerais) + 2:
            return

        try:
            novo_val = float(item.text().strip().replace(",", "."))
            if novo_val < 0 or not np.isfinite(novo_val):
                raise ValueError("Valor negativo ou não finito.")
        except ValueError:
            QMessageBox.warning(self, "Valor Inválido", "Insira apenas números reais não-negativos (ex: 0.045 ou 0,045).")
            self.carregar_tabela()
            return

        self.tabela.blockSignals(True)
        if col < len(minerais):
            mineral = minerais[col]
            KD_DATABASE[magma][mineral][el] = novo_val
            item.setText(f"{novo_val:.5f}")
        elif col == len(minerais):
            C0_DATABASE[el] = max(1e-12, novo_val)
            item.setText(f"{C0_DATABASE[el]:.4g}")
            if self.callback_c0_custom:
                self.callback_c0_custom()
        elif col == len(minerais) + 1:
            CA_PPM_DATABASE[el] = max(0.0, novo_val)
            item.setText(f"{CA_PPM_DATABASE[el]:.4g}")
            if self.callback_ca_custom:
                self.callback_ca_custom()

        razao_atual = CA_PPM_DATABASE[el] / max(1e-12, C0_DATABASE[el])
        item_raz = self.tabela.item(row, len(minerais) + 2)
        if item_raz:
            item_raz.setText(f"{razao_atual:.3f}")
        self.tabela.blockSignals(False)

        if self.callback_atualizar:
            self.callback_atualizar()

    def restaurar_padrao(self):
        global KD_DATABASE, C0_DATABASE, CA_PPM_DATABASE
        KD_DATABASE = copy.deepcopy(DEFAULT_KD_DATABASE)
        C0_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM['Manto Primitivo (Sun & McDonough, 1989)'])
        CA_PPM_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM['Crosta Continental Superior (Rudnick & Gao, 2003)'])
        self.carregar_tabela()
        if self.callback_reset:
            self.callback_reset()
        if self.callback_atualizar:
            self.callback_atualizar()
        QMessageBox.information(self, "Restaurado", "O banco de dados padrão foi restaurado com sucesso.")

    def exportar_csv(self):
        magma = self.cb_magma_tab.currentText()
        arquivo, _ = QFileDialog.getSaveFileName(
            self, f"Exportar Parâmetros ({magma})", f"parametros_{magma.split()[0].lower()}.csv", "CSV Files (*.csv)"
        )
        if arquivo:
            try:
                minerais = list(KD_DATABASE[magma].keys())
                dados = {"Elemento": ELEMENTOS_MASTER}
                for m in minerais:
                    dados[m] = [KD_DATABASE[magma][m][el] for el in ELEMENTOS_MASTER]
                dados["C0_ppm"] = [C0_DATABASE[el] for el in ELEMENTOS_MASTER]
                dados["Ca_ppm"] = [CA_PPM_DATABASE[el] for el in ELEMENTOS_MASTER]
                df = pd.DataFrame(dados)
                df.to_csv(arquivo, index=False, encoding='utf-8-sig')
                QMessageBox.information(self, "Sucesso", f"Tabela exportada para:\n{arquivo}")
            except Exception as e:
                QMessageBox.critical(self, "Erro", f"Falha ao exportar CSV:\n{str(e)}")

    @staticmethod
    def _converter_float_seguro(valor):
        if isinstance(valor, pd.Series):
            valor = valor.iloc[0]
        if pd.isna(valor):
            return None
        if isinstance(valor, (int, float, np.number)):
            val_f = float(valor)
            return val_f if np.isfinite(val_f) else None
        texto = str(valor).strip().replace(",", ".")
        if not texto:
            return None
        val_f = float(texto)
        return val_f if np.isfinite(val_f) else None

    @staticmethod
    def _detectar_delimitador_csv(caminho_arquivo):
        with open(caminho_arquivo, 'r', encoding='utf-8-sig', errors='replace') as f:
            for linha in f:
                linha_limpa = linha.strip()
                if linha_limpa and not linha_limpa.startswith('#'):
                    if ';' in linha_limpa:
                        return ';'
                    elif '\t' in linha_limpa:
                        return '\t'
                    return ','
        return ','

    def importar_csv(self):
        magma = self.cb_magma_tab.currentText()
        arquivo, _ = QFileDialog.getOpenFileName(
            self, f"Importar Parâmetros para {magma}", "", "CSV Files (*.csv);;All Files (*)"
        )
        if arquivo:
            try:
                delim = self._detectar_delimitador_csv(arquivo)
                df = pd.read_csv(arquivo, sep=delim, comment='#', encoding='utf-8-sig')
                df.columns = [str(c).strip() for c in df.columns]

                if "Elemento" not in df.columns:
                    raise ValueError("O arquivo CSV deve conter uma coluna chamada 'Elemento'.")

                df["Elemento"] = df["Elemento"].astype(str).str.strip()
                df = df.set_index("Elemento")
                df = df[~df.index.duplicated(keep='first')]

                minerais = list(KD_DATABASE[magma].keys())
                for m in minerais:
                    if m in df.columns:
                        for el in ELEMENTOS_MASTER:
                            if el in df.index:
                                val_conv = self._converter_float_seguro(df.loc[el, m])
                                if val_conv is not None and val_conv >= 0.0:
                                    KD_DATABASE[magma][m][el] = val_conv

                col_c0 = next((c for c in ("C0_ppm", "C0 (ppm)", "C0") if c in df.columns), None)
                if col_c0:
                    for el in ELEMENTOS_MASTER:
                        if el in df.index:
                            val_conv = self._converter_float_seguro(df.loc[el, col_c0])
                            if val_conv is not None and val_conv > 0.0:
                                C0_DATABASE[el] = val_conv
                    if self.callback_c0_custom:
                        self.callback_c0_custom()

                col_ca_ppm = next((c for c in ("Ca_ppm", "Ca (ppm)", "Ca") if c in df.columns), None)
                if col_ca_ppm:
                    for el in ELEMENTOS_MASTER:
                        if el in df.index:
                            val_conv = self._converter_float_seguro(df.loc[el, col_ca_ppm])
                            if val_conv is not None and val_conv >= 0.0:
                                CA_PPM_DATABASE[el] = val_conv
                    if self.callback_ca_custom:
                        self.callback_ca_custom()

                self.carregar_tabela()
                if self.callback_atualizar:
                    self.callback_atualizar()
                QMessageBox.information(self, "Sucesso", f"Dados importados com sucesso para o sistema {magma}.")
            except Exception as e:
                QMessageBox.critical(self, "Erro na Importação", f"Não foi possível processar o arquivo CSV:\n{str(e)}")


class TeoriaDialog(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Fundamentação Teórica e Equacionamento")
        self.resize(850, 730)
        layout = QVBoxLayout()
        texto_teoria = QTextBrowser()
        texto_teoria.setOpenExternalLinks(True)
        texto_teoria.setHtml("""
        <h2 style='color: #8B0000;'>Equacionamento Termodinâmico e Balanço de Massa</h2>
        <p>Nos modelos de fusão parcial (Shaw, 1970), define-se o coeficiente de partição inicial da fonte como
        <b>D<sub>0</sub> = &Sigma;(X<sub>0,i</sub> &times; Kd<sub>i</sub>)</b> e o parâmetro da reação de fusão como
        <b>P = &Sigma;(p<sub>i</sub> &times; Kd<sub>i</sub>)</b>, onde <i>X<sub>0,i</sub></i> é a fração modal em peso do mineral <i>i</i> na fonte
        (&Sigma;X<sub>0,i</sub> = 1) e <i>p<sub>i</sub></i> é a fração estequiométrica líquida com que o mineral <i>i</i> entra no fundido (&Sigma;p<sub>i</sub> = 1).
        A concentração absoluta no líquido (em ppm) é obtida multiplicando o fator de enriquecimento pela concentração inicial <i>C<sub>0</sub></i> (ppm).</p>

        <ul>
            <li><b>Fusão Modal (p<sub>i</sub> = X<sub>0,i</sub>):</b> Implica <b>P = D<sub>0</sub></b> e ausência de exaustão precoce de fases até F = 1.</li>
            <li><b>Fusão Não-Modal Cotética (p<sub>i</sub> &ge; 0, p<sub>i</sub> &ne; X<sub>0,i</sub>):</b> Todas as fases reagentes são consumidas em taxas distintas de suas modas iniciais.</li>
            <li><b>Fusão Não-Modal Peritética / Incongruente (p<sub>i</sub> &lt; 0 para fases peritéticas):</b> Uma fase mineral é precipitada enquanto as demais fundem (ex: fusão de lherzolito a espinélio gerando olivina peritética com <i>p<sub>Ol</sub> &lt; 0</i>, ou fusão por desidratação de biotita/anfibólio gerando granada/Opx peritético).</li>
            <li><b>Limite de Exaustão Mineral (F<sub>max</sub>):</b> A fração de massa residual de cada fase no sistema é <i>m<sub>i</sub>(F) = X<sub>0,i</sub> - p<sub>i</sub>F</i>. O estágio invariante/divariante encerra-se quando a primeira fase consumida (<i>p<sub>i</sub> &gt; 0</i>) atinge massa zero em <b>F<sub>max</sub> = min(X<sub>0,i</sub> / p<sub>i</sub>)</b>.</li>
            <li><b>Normalização ao Condrito C1 vs. Concentração Bruta (ppm):</b> A visualização direta em ppm bruto no diagrama multielementar exibe o efeito Oddo-Harkins (alternância em dente de serra entre elementos de número atômico par e ímpar). O modo <i>Normalizado ao Condrito C1 (Cl / C_Condrito)</i> divide a concentração absoluta calculada em ppm pelos valores do Condrito C1 (Sun & McDonough, 1989; McDonough & Sun, 1995), suavizando o efeito Oddo-Harkins e fornecendo as razões petrogenéticas reais <b>(La/Yb)<sub>N</sub></b> e <b>Eu/Eu*</b>.</li>
        </ul>

        <h3>1. Fusão em Equilíbrio Modal e Não-Modal (Batch Melting - Shaw, 1970)</h3>
        <p style='background-color: #f0f0f0; padding: 10px;'><b>C<sub>l</sub> / C<sub>0</sub> = 1 / [D<sub>0</sub> + F(1 - P)]</b> &nbsp;&rArr;&nbsp; <b>C<sub>l</sub> (ppm) = C<sub>0</sub> / [D<sub>0</sub> + F(1 - P)]</b></p>

        <h3>2. Fusão Fracionada Acumulada e Instantânea (Fractional Melting - Shaw, 1970)</h3>
        <p style='background-color: #f0f0f0; padding: 10px;'>
        <i>Líquido Acumulado (Pooled Melt):</i> <b>C&#772;<sub>l</sub> / C<sub>0</sub> = (1 / F) &times; [1 - (1 - P &times; F / D<sub>0</sub>)<sup>(1/P)</sup>]</b><br>
        <i>Líquido Instantâneo:</i> <b>C<sub>l</sub> / C<sub>0</sub> = (1 / D<sub>0</sub>) &times; (1 - P &times; F / D<sub>0</sub>)<sup>(1/P - 1)</sup></b><br>
        <i>Limite analítico quando P &rarr; 0:</i> <b>(1 - P &times; F / D<sub>0</sub>)<sup>(1/P)</sup> &rarr; exp(-F / D<sub>0</sub>)</b></p>

        <h3>3. Cristalização Fracionada (FC - Destilação de Rayleigh)</h3>
        <p style='background-color: #f0f0f0; padding: 10px;'><b>C<sub>l</sub> / C<sub>0</sub> = F<sup>(D - 1)</sup></b> &nbsp;&rArr;&nbsp; <b>C<sub>l</sub> (ppm) = C<sub>0</sub> &times; F<sup>(D - 1)</sup></b></p>

        <h3>4. Assimilação e Cristalização Fracionada (AFC - DePaolo, 1981)</h3>
        <p>Sendo <i>r = M&#775;<sub>a</sub> / M&#775;<sub>c</sub></i> (com <i>r &lt; 1</i>) e <i>z = (r + D - 1) / (r - 1)</i>:</p>
        <p style='background-color: #f0f0f0; padding: 10px;'><b>C<sub>l</sub> / C<sub>0</sub> = F<sup>-z</sup> + [r / (r + D - 1)] &times; (C<sub>a</sub> / C<sub>0</sub>) &times; (1 - F<sup>-z</sup>)</b><br>
        <i>Em unidades absolutas (ppm):</i> <b>C<sub>l</sub> = C<sub>0</sub> F<sup>-z</sup> + [r / (r + D - 1)] &times; C<sub>a</sub> &times; (1 - F<sup>-z</sup>)</b><br>
        <i>Limite analítico em z &rarr; 0 (r + D = 1):</i> <b>C<sub>l</sub> / C<sub>0</sub> = 1 + [r / (r - 1)] &times; (C<sub>a</sub> / C<sub>0</sub>) &times; ln(F)</b><br>
        <i>Balanço de Massa Assimilada Acumulada:</i> <b>M<sub>a</sub> / M<sub>m</sub><sup>0</sup> = [r / (1 - r)] &times; (1 - F)</b></p>
        """)
        layout.addWidget(texto_teoria)
        btn_fechar = QPushButton("Fechar")
        btn_fechar.clicked.connect(self.accept)
        layout.addWidget(btn_fechar, alignment=Qt.AlignRight)
        self.setLayout(layout)


class TelaBoasVindas(QDialog):
    def __init__(self, pixmap_icone=None):
        super().__init__()
        self.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        self.setFixedSize(680, 410)
        self.setStyleSheet("QDialog { background-color: #181818; border: 2px solid #8B0000; }")

        layout = QVBoxLayout()
        layout.setAlignment(Qt.AlignCenter)

        if pixmap_icone is not None:
            lbl_ico = QLabel()
            lbl_ico.setPixmap(pixmap_icone.scaled(76, 76, Qt.KeepAspectRatio, Qt.SmoothTransformation))
            lbl_ico.setAlignment(Qt.AlignCenter)
            layout.addWidget(lbl_ico)

        lbl_logo = QLabel("Magma_Trace")
        lbl_logo.setStyleSheet("font-family: 'Segoe UI', Arial; font-size: 48px; font-weight: bold; color: #e64a19;")
        lbl_logo.setAlignment(Qt.AlignCenter)

        lbl_sub = QLabel("Modelagem Quantitativa de Elementos-Traço em Processos Ígneos")
        lbl_sub.setStyleSheet("font-size: 16px; color: #cccccc; margin-bottom: 16px;")
        lbl_sub.setAlignment(Qt.AlignCenter)

        lbl_creditos = QLabel(
            "<b>Magma_Trace (C) 2026</b><br>"
            "Prof. Dr. Felipe P. Leitzke<br>"
            "Departamento de Mineralogia e Petrologia — IGEO / UFRGS"
        )
        lbl_creditos.setStyleSheet("font-size: 14px; color: #9e9e9e; margin-top: 10px;")
        lbl_creditos.setAlignment(Qt.AlignCenter)

        btn_iniciar = QPushButton("Iniciar Programa")
        btn_iniciar.setCursor(Qt.PointingHandCursor)
        btn_iniciar.setStyleSheet("""
            QPushButton { background-color: #8B0000; color: white; font-size: 15px; font-weight: bold; border-radius: 3px; padding: 10px 28px; margin-top: 25px; }
            QPushButton:hover { background-color: #b71c1c; }
        """)
        btn_iniciar.clicked.connect(self.accept)

        layout.addWidget(lbl_logo)
        layout.addWidget(lbl_sub)
        layout.addWidget(lbl_creditos)
        layout.addWidget(btn_iniciar, alignment=Qt.AlignCenter)
        self.setLayout(layout)


class GeochemApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Magma_Trace — Modelagem Geoquímica de Elementos-Traço")
        self.resize(1380, 860)

        self.em_ppm_plot = False
        self.norm_condrito_plot = False
        self.modo_spider_plot = False
        self.F_plot = None
        self.resultados_plot = None
        self.F_spider_plot = None
        self.res_spider_plot = None
        self.el_ordenados_spider = []
        self.elementos_plot = []
        self.processo_plot = ""
        self.magma_plot = ""
        self.bulk_D0_plot = {}
        self.bulk_P_plot = {}
        self.moda_x0_plot = {}
        self.moda_p_plot = {}
        self.f_max_plot = 1.0
        self.min_exauridos_plot = []

        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        layout_principal = QHBoxLayout(main_widget)

        painel_controles = QVBoxLayout()
        painel_controles.setContentsMargins(5, 5, 12, 5)

        grupo_config = QGroupBox("Configuração do Sistema e Visualização")
        form_config = QFormLayout()

        self.cb_processo = QComboBox()
        self.cb_processo.addItems([
            "Fusão em Equilíbrio (Batch Melting - Shaw, 1970)",
            "Fusão Fracionada Acumulada (Pooled Fractional Melting)",
            "Fusão Fracionada Instantânea (Instantaneous Melt)",
            "Cristalização Fracionada (FC - Rayleigh)",
            "Assimilação + Fracionamento (AFC - DePaolo, 1981)"
        ])
        self.cb_processo.currentIndexChanged.connect(self.atualizar_ui)

        self.cb_magma = QComboBox()
        self.cb_magma.addItems(list(KD_DATABASE.keys()))
        self.cb_magma.currentIndexChanged.connect(self.ao_mudar_sistema_magmatico)

        self.cb_modo_grafico = QComboBox()
        self.cb_modo_grafico.addItems([
            "1. Evolução Relativa Contínua (Cl/C0 vs. F)",
            "2. Diagrama Multielementar Relativo (Cl/C0 vs. Elemento)",
            "3. Concentração Absoluta Contínua (Cl em ppm vs. F)",
            "4. Diagrama Multielementar Normalizado ao Condrito C1 (Cl/C_Condrito)",
            "5. Diagrama Multielementar em Concentração Bruta (Cl em ppm)"
        ])
        self.cb_modo_grafico.currentIndexChanged.connect(self.atualizar_ui)

        self.cb_preset_c0 = QComboBox()
        self.cb_preset_c0.addItems(list(RESERVATORIOS_C0_PPM.keys()) + ["Personalizado (Definido na Tabela)"])
        self.cb_preset_c0.currentIndexChanged.connect(self.selecionar_reservatorio_c0)
        self.lbl_preset_c0 = QLabel("Reservatório C0:")

        self.chk_clamp_y = QCheckBox("Restringir envelope logarítmico do Eixo Y")
        self.chk_clamp_y.setChecked(True)
        self.chk_clamp_y.setToolTip("Impede que underflow numérico em elementos extremos comprima o eixo logarítmico.")
        self.chk_clamp_y.stateChanged.connect(self.atualizar_grafico)

        form_config.addRow("Processo:", self.cb_processo)
        form_config.addRow("Sistema:", self.cb_magma)
        form_config.addRow("Gráfico:", self.cb_modo_grafico)
        form_config.addRow(self.lbl_preset_c0, self.cb_preset_c0)
        form_config.addRow("", self.chk_clamp_y)
        grupo_config.setLayout(form_config)
        painel_controles.addWidget(grupo_config)

        grupo_elementos = QGroupBox("Elementos-Traço Ativos")
        layout_elementos = QVBoxLayout()

        self.lista_elementos = QListWidget()
        self.lista_elementos.setFixedHeight(105)
        elementos_padrao = ['Rb', 'Sr', 'Nd', 'Sm', 'Eu', 'Yb']

        for el in ELEMENTOS_MASTER:
            item = QListWidgetItem(el)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if el in elementos_padrao else Qt.Unchecked)
            self.lista_elementos.addItem(item)

        self.lista_elementos.itemChanged.connect(self.atualizar_grafico)
        layout_elementos.addWidget(self.lista_elementos)

        botoes_sel_el = QHBoxLayout()
        btn_sel_etr = QPushButton("ETR (La-Lu)")
        btn_sel_etr.setToolTip("Selecionar os 14 Elementos Terras Raras")
        btn_sel_etr.clicked.connect(lambda: self.selecionar_conjunto_elementos(LISTA_ETR))

        btn_sel_spider = QPushButton("Spidergram")
        btn_sel_spider.setToolTip("Selecionar LILE + HFSE + ETR")
        btn_sel_spider.clicked.connect(lambda: self.selecionar_conjunto_elementos(ORDEM_SPIDERGRAM[:-4]))

        btn_sel_padrao = QPushButton("Padrão")
        btn_sel_padrao.clicked.connect(lambda: self.selecionar_conjunto_elementos(['Rb', 'Sr', 'Nd', 'Sm', 'Eu', 'Yb']))

        btn_sel_limpar = QPushButton("Limpar")
        btn_sel_limpar.clicked.connect(lambda: self.selecionar_conjunto_elementos([]))

        for b in (btn_sel_etr, btn_sel_spider, btn_sel_padrao, btn_sel_limpar):
            b.setStyleSheet("padding: 3px 6px; font-size: 11px;")
            botoes_sel_el.addWidget(b)

        layout_elementos.addLayout(botoes_sel_el)
        grupo_elementos.setLayout(layout_elementos)
        painel_controles.addWidget(grupo_elementos)

        self.grupo_minerais = QGroupBox("Proporções Mineralógicas: Fonte (X0) e Reação de Fusão (pi)")
        layout_min_container = QVBoxLayout()

        self.chk_peritetico = QCheckBox("Habilitar Fusão Peritética / Incongruente (permitir pi < 0)")
        self.chk_peritetico.setChecked(False)
        self.chk_peritetico.stateChanged.connect(self.alternar_modo_peritetico)
        layout_min_container.addWidget(self.chk_peritetico)

        self.grid_minerais = QGridLayout()
        self.grid_minerais.setHorizontalSpacing(8)
        self.grid_minerais.setVerticalSpacing(3)

        self.lbl_col_min = QLabel("<b>Fase</b>")
        self.lbl_col_x0 = QLabel("<b>Moda X0 (%)</b>")
        self.lbl_col_x0.setAlignment(Qt.AlignCenter)
        self.lbl_col_p = QLabel("<b>Reação pi (%)</b>")
        self.lbl_col_p.setAlignment(Qt.AlignCenter)

        self.grid_minerais.addWidget(self.lbl_col_min, 0, 0)
        self.grid_minerais.addWidget(self.lbl_col_x0, 0, 1, 1, 2)
        self.grid_minerais.addWidget(self.lbl_col_p, 0, 3, 1, 2)

        self.sliders_x0 = {}
        self.labels_x0 = {}
        self.sliders_p = {}
        self.labels_p = {}
        self.row_labels = {}

        moda_inicial = MODA_PADRAO_SISTEMA['Basáltico (Máfico)']

        for idx, min_nome in enumerate(TODOS_MINERAIS, start=1):
            lbl_nome = QLabel(min_nome + ":")
            self.row_labels[min_nome] = lbl_nome
            val_ini = moda_inicial.get(min_nome, 0)

            sl_x = QSlider(Qt.Horizontal)
            sl_x.setRange(0, 100)
            sl_x.setValue(val_ini)
            sl_x.valueChanged.connect(self.atualizar_grafico)
            lbl_x = QLabel(f"{val_ini}%")
            lbl_x.setFixedWidth(36)
            lbl_x.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

            sl_p = QSlider(Qt.Horizontal)
            sl_p.setRange(0, 100)
            sl_p.setValue(val_ini)
            sl_p.valueChanged.connect(self.atualizar_grafico)
            lbl_p = QLabel(f"{val_ini}%")
            lbl_p.setFixedWidth(38)
            lbl_p.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

            self.sliders_x0[min_nome] = sl_x
            self.labels_x0[min_nome] = lbl_x
            self.sliders_p[min_nome] = sl_p
            self.labels_p[min_nome] = lbl_p

            self.grid_minerais.addWidget(lbl_nome, idx, 0)
            self.grid_minerais.addWidget(sl_x, idx, 1)
            self.grid_minerais.addWidget(lbl_x, idx, 2)
            self.grid_minerais.addWidget(sl_p, idx, 3)
            self.grid_minerais.addWidget(lbl_p, idx, 4)

        layout_min_container.addLayout(self.grid_minerais)

        self.btn_copiar_moda = QPushButton("Igualar pi = X0 (Fusão Modal: P = D0)")
        self.btn_copiar_moda.clicked.connect(self.igualar_p_a_x0)
        layout_min_container.addWidget(self.btn_copiar_moda)

        self.lbl_soma_x0 = QLabel("<b>Soma X0: 100%</b>")
        self.lbl_soma_p = QLabel("<b>Soma pi: 100%</b>")
        layout_somas = QHBoxLayout()
        layout_somas.addWidget(self.lbl_soma_x0)
        layout_somas.addWidget(self.lbl_soma_p)
        layout_min_container.addLayout(layout_somas)

        self.grupo_minerais.setLayout(layout_min_container)
        painel_controles.addWidget(self.grupo_minerais)

        self.grupo_afc = QGroupBox("Parâmetros de Assimilação (AFC - DePaolo, 1981)")
        layout_afc = QFormLayout()

        self.sl_r = QSlider(Qt.Horizontal)
        self.sl_r.setRange(1, 95)
        self.sl_r.setValue(20)
        self.sl_r.valueChanged.connect(self.atualizar_grafico)
        self.lbl_r = QLabel("0.20")
        self.lbl_r.setFixedWidth(40)
        row_r = QHBoxLayout()
        row_r.addWidget(self.sl_r)
        row_r.addWidget(self.lbl_r)
        layout_afc.addRow("Taxa r (Ma/Mc):", row_r)

        self.cb_modo_ca = QComboBox()
        self.cb_modo_ca.addItems([
            "Por Elemento (Razão Real Ca_ppm / C0_ppm)",
            "Paramétrico Global (Razão Ca/C0 fixa no Slider)"
        ])
        self.cb_modo_ca.currentIndexChanged.connect(self.atualizar_ui)
        layout_afc.addRow("Modo Assimilante:", self.cb_modo_ca)

        self.cb_preset_ca = QComboBox()
        self.cb_preset_ca.addItems(list(RESERVATORIOS_C0_PPM.keys()) + ["Personalizado (Definido na Tabela)"])
        self.cb_preset_ca.setCurrentIndex(3)  # Crosta Continental Superior (Rudnick & Gao, 2003)
        self.cb_preset_ca.currentIndexChanged.connect(self.selecionar_reservatorio_ca)
        self.lbl_preset_ca = QLabel("Reservatório Ca:")
        layout_afc.addRow(self.lbl_preset_ca, self.cb_preset_ca)

        self.sl_ca_c0 = QSlider(Qt.Horizontal)
        self.sl_ca_c0.setRange(1, 100)
        self.sl_ca_c0.setValue(20)
        self.sl_ca_c0.valueChanged.connect(self.atualizar_grafico)
        self.lbl_ca_c0 = QLabel("2.00")
        self.lbl_ca_c0.setFixedWidth(40)
        self.lbl_titulo_ca_slider = QLabel("Ca/C0 (Global):")
        row_ca = QHBoxLayout()
        row_ca.addWidget(self.sl_ca_c0)
        row_ca.addWidget(self.lbl_ca_c0)
        layout_afc.addRow(self.lbl_titulo_ca_slider, row_ca)

        self.lbl_balanco_afc = QLabel("")
        self.lbl_balanco_afc.setWordWrap(True)
        layout_afc.addRow(self.lbl_balanco_afc)

        self.grupo_afc.setLayout(layout_afc)
        painel_controles.addWidget(self.grupo_afc)

        painel_controles.addStretch()

        botoes_layout1 = QHBoxLayout()
        btn_tabela = QPushButton("Parâmetros (Kd / C0 / Ca)")
        btn_tabela.clicked.connect(self.abrir_tabela_kd)
        btn_exportar = QPushButton("Exportar CSV")
        btn_exportar.clicked.connect(self.exportar_dados)
        btn_salvar_fig = QPushButton("Salvar Figura")
        btn_salvar_fig.clicked.connect(self.salvar_figura)
        botoes_layout1.addWidget(btn_tabela)
        botoes_layout1.addWidget(btn_exportar)
        botoes_layout1.addWidget(btn_salvar_fig)
        painel_controles.addLayout(botoes_layout1)

        botoes_layout2 = QHBoxLayout()
        btn_teoria = QPushButton("Equacionamento")
        btn_teoria.clicked.connect(self.abrir_teoria)
        btn_refs = QPushButton("Bibliografia")
        btn_refs.clicked.connect(self.mostrar_referencias)
        botoes_layout2.addWidget(btn_teoria)
        botoes_layout2.addWidget(btn_refs)
        painel_controles.addLayout(botoes_layout2)

        lbl_creditos = QLabel("<small><b>Magma_Trace (C) 2026</b> — Felipe P. Leitzke (DEMIPE/IGEO/UFRGS)</small>")
        lbl_creditos.setAlignment(Qt.AlignCenter)
        lbl_creditos.setStyleSheet("margin-top: 4px; color: #555555;")
        painel_controles.addWidget(lbl_creditos)

        self.figura = Figure(figsize=(8.5, 6.2), dpi=100)
        self.figura.subplots_adjust(left=0.10, right=0.96, top=0.91, bottom=0.11)
        self.canvas = FigureCanvas(self.figura)
        self.ax = self.figura.add_subplot(111)

        layout_principal.addLayout(painel_controles, 2)
        layout_principal.addWidget(self.canvas, 5)

        self.atualizar_ui()

    def ao_mudar_sistema_magmatico(self):
        magma_atual = self.cb_magma.currentText()
        padrao_sistema = MODA_PADRAO_SISTEMA.get(magma_atual, {})
        for min_nome in TODOS_MINERAIS:
            v = padrao_sistema.get(min_nome, 0)
            self.sliders_x0[min_nome].blockSignals(True)
            self.sliders_p[min_nome].blockSignals(True)
            self.sliders_x0[min_nome].setValue(v)
            self.sliders_p[min_nome].setValue(v)
            self.labels_x0[min_nome].setText(f"{v}%")
            self.labels_p[min_nome].setText(f"{v}%")
            self.sliders_x0[min_nome].blockSignals(False)
            self.sliders_p[min_nome].blockSignals(False)
        self.atualizar_ui()

    def selecionar_reservatorio_c0(self):
        global C0_DATABASE
        nome_res = self.cb_preset_c0.currentText()
        if nome_res in RESERVATORIOS_C0_PPM:
            C0_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM[nome_res])
            self.atualizar_grafico()

    def selecionar_reservatorio_ca(self):
        global CA_PPM_DATABASE
        nome_res = self.cb_preset_ca.currentText()
        if nome_res in RESERVATORIOS_C0_PPM:
            CA_PPM_DATABASE = copy.deepcopy(RESERVATORIOS_C0_PPM[nome_res])
            self.atualizar_grafico()

    def marcar_c0_personalizado(self):
        idx_custom = self.cb_preset_c0.count() - 1
        self.cb_preset_c0.blockSignals(True)
        self.cb_preset_c0.setCurrentIndex(idx_custom)
        self.cb_preset_c0.blockSignals(False)

    def marcar_ca_personalizado(self):
        idx_custom = self.cb_preset_ca.count() - 1
        self.cb_preset_ca.blockSignals(True)
        self.cb_preset_ca.setCurrentIndex(idx_custom)
        self.cb_preset_ca.blockSignals(False)

    def resetar_combos_padrao(self):
        self.cb_preset_c0.blockSignals(True)
        self.cb_preset_c0.setCurrentIndex(0)
        self.cb_preset_c0.blockSignals(False)

        self.cb_preset_ca.blockSignals(True)
        self.cb_preset_ca.setCurrentIndex(3)
        self.cb_preset_ca.blockSignals(False)

    def selecionar_conjunto_elementos(self, conjunto):
        self.lista_elementos.blockSignals(True)
        for i in range(self.lista_elementos.count()):
            item = self.lista_elementos.item(i)
            item.setCheckState(Qt.Checked if item.text() in conjunto else Qt.Unchecked)
        self.lista_elementos.blockSignals(False)
        self.atualizar_grafico()

    def limpar_memoria_plot(self):
        self.modo_spider_plot = False
        self.em_ppm_plot = False
        self.norm_condrito_plot = False
        self.F_plot = None
        self.resultados_plot = None
        self.F_spider_plot = None
        self.res_spider_plot = None
        self.el_ordenados_spider = []
        self.elementos_plot = []
        self.processo_plot = ""
        self.magma_plot = ""
        self.bulk_D0_plot = {}
        self.bulk_P_plot = {}
        self.moda_x0_plot = {}
        self.moda_p_plot = {}
        self.f_max_plot = 1.0
        self.min_exauridos_plot = []

    def alternar_modo_peritetico(self):
        permitir_neg = self.chk_peritetico.isChecked()
        min_val = -100 if permitir_neg else 0
        for min_nome in TODOS_MINERAIS:
            sl_p = self.sliders_p[min_nome]
            sl_p.blockSignals(True)
            if not permitir_neg and sl_p.value() < 0:
                sl_p.setValue(0)
                self.labels_p[min_nome].setText("0%")
            sl_p.setRange(min_val, 100)
            sl_p.blockSignals(False)
        self.atualizar_grafico()

    def igualar_p_a_x0(self):
        magma_atual = self.cb_magma.currentText()
        for min_nome in KD_DATABASE[magma_atual].keys():
            val_x0 = self.sliders_x0[min_nome].value()
            self.sliders_p[min_nome].blockSignals(True)
            self.sliders_p[min_nome].setValue(val_x0)
            self.labels_p[min_nome].setText(f"{val_x0}%")
            self.sliders_p[min_nome].blockSignals(False)
        self.atualizar_grafico()

    def atualizar_ui(self):
        magma_atual = self.cb_magma.currentText()
        processo = self.cb_processo.currentText()
        texto_modo = self.cb_modo_grafico.currentText()
        minerais_permitidos = KD_DATABASE[magma_atual].keys()

        is_fusao = "Fusão" in processo
        is_afc = "AFC" in processo
        usar_global_afc = "Global" in self.cb_modo_ca.currentText()
        usa_c0_ppm = ("ppm" in texto_modo) or ("Condrito" in texto_modo) or (is_afc and not usar_global_afc)

        self.lbl_preset_c0.setVisible(usa_c0_ppm)
        self.cb_preset_c0.setVisible(usa_c0_ppm)

        self.lbl_col_x0.setText("<b>Moda X0 (%)</b>" if is_fusao else "<b>Assembleia Fracionante (%)</b>")
        self.chk_peritetico.setVisible(is_fusao)
        self.lbl_col_p.setVisible(is_fusao)
        self.btn_copiar_moda.setVisible(is_fusao)
        self.lbl_soma_p.setVisible(is_fusao)

        for min_nome in TODOS_MINERAIS:
            visivel = min_nome in minerais_permitidos
            self.row_labels[min_nome].setVisible(visivel)
            self.sliders_x0[min_nome].setVisible(visivel)
            self.labels_x0[min_nome].setVisible(visivel)
            self.sliders_p[min_nome].setVisible(visivel and is_fusao)
            self.labels_p[min_nome].setVisible(visivel and is_fusao)

            if not visivel:
                self.sliders_x0[min_nome].blockSignals(True)
                self.sliders_x0[min_nome].setValue(0)
                self.labels_x0[min_nome].setText("0%")
                self.sliders_x0[min_nome].blockSignals(False)

                self.sliders_p[min_nome].blockSignals(True)
                self.sliders_p[min_nome].setValue(0)
                self.labels_p[min_nome].setText("0%")
                self.sliders_p[min_nome].blockSignals(False)

        self.grupo_afc.setVisible(is_afc)
        if is_afc:
            self.lbl_preset_ca.setVisible(not usar_global_afc)
            self.cb_preset_ca.setVisible(not usar_global_afc)
            self.lbl_titulo_ca_slider.setVisible(usar_global_afc)
            self.sl_ca_c0.setVisible(usar_global_afc)
            self.lbl_ca_c0.setVisible(usar_global_afc)

        self.atualizar_grafico()

    def abrir_tabela_kd(self):
        dialog = TabelaKdDialog(
            callback_atualizar=self.atualizar_grafico,
            callback_c0_custom=self.marcar_c0_personalizado,
            callback_ca_custom=self.marcar_ca_personalizado,
            callback_reset=self.resetar_combos_padrao
        )
        dialog.exec_()

    def abrir_teoria(self):
        dialog = TeoriaDialog()
        dialog.exec_()

    def mostrar_referencias(self):
        texto = (
            "<h3>Bibliografia Fundamental</h3>"
            "<p>- <i>Arth, J. G. (1976). Behavior of trace elements during magmatic processes: a summary of theoretical models and their applications. Journal of Research of the U.S. Geological Survey, 4(1), 41-47.</i><br>"
            "- <i>DePaolo, D. J. (1981). Trace element and isotopic effects of combined wallrock assimilation and fractional crystallization. Earth and Planetary Science Letters, 53(2), 189-202.</i><br>"
            "- <i>Fujimaki, H., Tatsumoto, M., & Aoki, K. I. (1984). Partition coefficients of Hf, Zr, and REE between phenocrysts and groundmasses. Journal of Geophysical Research: Solid Earth, 89(S02), B662-B672.</i><br>"
            "- <i>McDonough, W. F., & Sun, S. S. (1995). The composition of the Earth. Chemical Geology, 120(3-4), 223-253.</i><br>"
            "- <i>Rudnick, R. L., & Gao, S. (2003). Composition of the continental crust. Treatise on Geochemistry, 3, 1-64.</i><br>"
            "- <i>Shaw, D. M. (1970). Trace element fractionation during anatexis. Geochimica et Cosmochimica Acta, 34(2), 237-243.</i><br>"
            "- <i>Sun, S. S., & McDonough, W. F. (1989). Chemical and isotopic systematics of oceanic basalts: implications for mantle composition and processes. Geological Society, London, Special Publications, 42(1), 313-345.</i><br>"
            "- <i>Rollinson, H. R. (1993). Using Geochemical Data: Evaluation, Presentation, Interpretation. Longman Scientific & Technical.</i><br>"
            "- <i>Rollinson, H., & Pease, V. (2021). Using Geochemical Data: To Understand Geological Processes (2nd ed.). Cambridge University Press.</i></p>"
        )
        msg = QMessageBox(self)
        msg.setWindowTitle("Bibliografia")
        msg.setText(texto)
        msg.setIcon(QMessageBox.Information)
        msg.exec_()

    def salvar_figura(self):
        if not self.elementos_plot or self.resultados_plot is None:
            QMessageBox.warning(self, "Aviso", "Não há gráfico válido na tela para salvar.")
            return

        arquivo, _ = QFileDialog.getSaveFileName(
            self, "Salvar Figura de Publicação", "grafico_magma_trace.png",
            "PNG Image 300 DPI (*.png);;PDF Vector (*.pdf);;SVG Vector (*.svg)"
        )
        if arquivo:
            try:
                self.figura.savefig(arquivo, dpi=300, bbox_inches='tight')
                QMessageBox.information(self, "Sucesso", f"Figura salva em alta resolução em:\n{arquivo}")
            except Exception as e:
                QMessageBox.critical(self, "Erro", f"Falha ao salvar a figura:\n{str(e)}")

    def exportar_dados(self):
        if not self.elementos_plot or self.resultados_plot is None:
            QMessageBox.warning(self, "Aviso", "Não há dados válidos calculados para exportar.")
            return

        if self.norm_condrito_plot:
            unid_str = "Condrito_C1"
        elif self.em_ppm_plot:
            unid_str = "ppm"
        else:
            unid_str = "Cl_C0"

        prefixo = "spidergram" if self.modo_spider_plot else "evolucao_F"
        nome_sugerido = f"{prefixo}_{unid_str}_magma_trace.csv"

        opcoes = QFileDialog.Options()
        arquivo, _ = QFileDialog.getSaveFileName(
            self, "Exportar Resultados Magma_Trace", nome_sugerido,
            "CSV Files (*.csv);;Text Files (*.txt)", options=opcoes
        )

        if arquivo:
            try:
                is_fusao = "Fusão" in self.processo_plot
                is_afc = "AFC" in self.processo_plot
                usar_global = "Global" in self.cb_modo_ca.currentText()
                modo_str = self.cb_modo_grafico.currentText()
                usa_c0 = self.em_ppm_plot or self.norm_condrito_plot or (is_afc and not usar_global)

                if self.modo_spider_plot and self.res_spider_plot is not None:
                    dados_spider = {
                        "Elemento": self.el_ordenados_spider,
                        "Bulk_D0": [self.bulk_D0_plot[el] for el in self.el_ordenados_spider]
                    }
                    if is_fusao:
                        dados_spider["Parametro_P"] = [self.bulk_P_plot[el] for el in self.el_ordenados_spider]
                    if usa_c0:
                        dados_spider["C0_ppm"] = [C0_DATABASE[el] for el in self.el_ordenados_spider]
                    if self.norm_condrito_plot:
                        dados_spider["Condrito_C1_ppm"] = [CONDRITO_C1_PPM[el] for el in self.el_ordenados_spider]
                    if is_afc:
                        ca_val_g = self.sl_ca_c0.value() / 10.0
                        if not usar_global:
                            dados_spider["Ca_ppm"] = [CA_PPM_DATABASE[el] for el in self.el_ordenados_spider]
                        dados_spider["Ca_C0"] = [
                            ca_val_g if usar_global else (CA_PPM_DATABASE[el] / max(1e-12, C0_DATABASE[el]))
                            for el in self.el_ordenados_spider
                        ]

                    if self.norm_condrito_plot:
                        rotulo_col = "Cl_CondritoC1"
                    elif self.em_ppm_plot:
                        rotulo_col = "Cl_ppm"
                    else:
                        rotulo_col = "Cl_C0"

                    for idx_f, f_val in enumerate(self.F_spider_plot):
                        dados_spider[f"{rotulo_col} (F={f_val:.4g})"] = [self.res_spider_plot[el][idx_f] for el in self.el_ordenados_spider]

                    df = pd.DataFrame(dados_spider)
                else:
                    col_f = "F_fusao (Fracao_Liquido_Gerado)" if is_fusao else "F_liq (Fracao_Liquido_Remanescente)"
                    df = pd.DataFrame({col_f: self.F_plot})
                    sufixo_col = "(ppm)" if self.em_ppm_plot else "(Cl/C0)"
                    for el in self.elementos_plot:
                        df[f"{el} {sufixo_col}"] = self.resultados_plot[el]

                x0_str = ", ".join([f"{m}={p*100:.1f}%" for m, p in self.moda_x0_plot.items() if abs(p) > 0])
                d0_str = ", ".join([f"{el}={self.bulk_D0_plot[el]:.4f}" for el in self.elementos_plot])

                with open(arquivo, 'w', encoding='utf-8-sig', newline='') as f:
                    f.write("# --- MAGMA_TRACE: RELATORIO DE EXPORTACAO ---\n")
                    f.write(f"# Modo de Visualizacao Exportado: {modo_str}\n")
                    f.write(f"# Processo Modelado: {self.processo_plot}\n")
                    f.write(f"# Sistema Magmatico: {self.magma_plot}\n")
                    f.write(f"# Moda Normalizada (X0): {x0_str}\n")
                    f.write(f"# Bulk D0: {d0_str}\n")

                    if usa_c0:
                        c0_str = ", ".join([f"{el}={C0_DATABASE[el]:.4g}" for el in self.elementos_plot])
                        f.write(f"# Reservatorio C0: {self.cb_preset_c0.currentText()}\n")
                        f.write(f"# Concentracao Inicial C0 (ppm): {c0_str}\n")

                    if is_fusao:
                        p_str = ", ".join([f"{m}={p*100:.1f}%" for m, p in self.moda_p_plot.items() if abs(p) > 0])
                        P_el_str = ", ".join([f"{el}={self.bulk_P_plot[el]:.4f}" for el in self.elementos_plot])
                        f.write(f"# Proporcao de Reacao Normalizada (p_i): {p_str}\n")
                        f.write(f"# Parametro P (Shaw, 1970): {P_el_str}\n")
                        if self.min_exauridos_plot:
                            nomes_ex = ", ".join(self.min_exauridos_plot)
                            f.write(f"# Limite de Exaustao Mineral (F_max): {self.f_max_plot:.4f} (Exaustao de {nomes_ex})\n")

                    if is_afc:
                        r_val = self.sl_r.value() / 100.0
                        f.write(f"# Taxa r (Ma/Mc): {r_val:.2f}\n")
                        if usar_global:
                            f.write(f"# Razao Ca/C0 (Global): {self.sl_ca_c0.value() / 10.0:.2f}\n")
                        else:
                            ca_ppm_str = ", ".join([f"{el}={CA_PPM_DATABASE[el]:.4g}" for el in self.elementos_plot])
                            ca_rat_str = ", ".join([f"{el}={CA_PPM_DATABASE[el]/max(1e-12, C0_DATABASE[el]):.3f}" for el in self.elementos_plot])
                            f.write(f"# Reservatorio Assimilante Ca: {self.cb_preset_ca.currentText()}\n")
                            f.write(f"# Concentracao Assimilante Ca (ppm): {ca_ppm_str}\n")
                            f.write(f"# Razao Ca/C0 Resultante: {ca_rat_str}\n")
                    f.write("#\n")
                    df.to_csv(f, index=False)

                QMessageBox.information(self, "Sucesso", f"Matriz exportada com metadados para:\n{arquivo}")
            except Exception as e:
                QMessageBox.critical(self, "Erro na Exportação", f"Ocorreu uma falha ao salvar o arquivo:\n{str(e)}")

    def calcular_concentracoes_em_F(self, F_vec, elementos_ativos, bulk_D0, bulk_P, processo, r_afc, ca_c0_global, usar_ca_global, em_ppm=False, norm_condrito=False):
        resultados = {el: np.zeros_like(F_vec, dtype=float) for el in elementos_ativos}
        log_floor = np.log(1e-30)
        log_ceil = np.log(1e30)

        for el in elementos_ativos:
            D0 = bulk_D0[el]
            P = bulk_P[el]
            D0_eff = max(D0, 1e-12)
            c0_el = max(1e-12, C0_DATABASE.get(el, 1.0))
            ca_c0_el = ca_c0_global if usar_ca_global else (CA_PPM_DATABASE.get(el, 1.0) / c0_el)

            if norm_condrito:
                fator_escala = c0_el / CONDRITO_C1_PPM.get(el, 1.0)
            elif em_ppm:
                fator_escala = c0_el
            else:
                fator_escala = 1.0

            if "Batch Melting" in processo:
                denom = np.maximum(1e-15, D0_eff + F_vec * (1.0 - P))
                cl_c0 = np.clip(1.0 / denom, 1e-30, 1e30)

            elif "Fusão Fracionada Acumulada" in processo:
                if abs(P) < 1e-8:
                    res = -np.expm1(-F_vec / D0_eff) / F_vec
                else:
                    u = -(P * F_vec) / D0_eff
                    u_clipped = np.maximum(-1.0 + 1e-15, u)
                    numerador = -np.expm1((1.0 / P) * np.log1p(u_clipped))
                    res = numerador / F_vec
                cl_c0 = np.clip(res, 1e-30, 1e30)

            elif "Fusão Fracionada Instantânea" in processo:
                if abs(P) < 1e-8:
                    log_cl = -np.log(D0_eff) - (F_vec / D0_eff)
                else:
                    u = -(P * F_vec) / D0_eff
                    u_clipped = np.maximum(-1.0 + 1e-15, u)
                    log_cl = -np.log(D0_eff) + (1.0 / P - 1.0) * np.log1p(u_clipped)
                cl_c0 = np.exp(np.clip(log_cl, log_floor, log_ceil))

            elif "FC" in processo and "AFC" not in processo:
                log_cl = (D0 - 1.0) * np.log(F_vec)
                cl_c0 = np.exp(np.clip(log_cl, log_floor, log_ceil))

            elif "AFC" in processo:
                if abs(r_afc + D0 - 1.0) < 1e-8:
                    res = 1.0 + (r_afc / (r_afc - 1.0)) * ca_c0_el * np.log(F_vec)
                else:
                    z = (r_afc + D0 - 1.0) / (r_afc - 1.0)
                    exp_arg = np.clip(-z * np.log(F_vec), log_floor, log_ceil)
                    termo_f = np.exp(exp_arg)
                    termo_assim = -(r_afc / (r_afc + D0 - 1.0)) * ca_c0_el * np.expm1(exp_arg)
                    res = termo_f + termo_assim
                cl_c0 = np.clip(res, 1e-30, 1e30)

            resultados[el] = cl_c0 * fator_escala
        return resultados

    def atualizar_grafico(self):
        magma_atual = self.cb_magma.currentText()
        processo = self.cb_processo.currentText()
        texto_modo = self.cb_modo_grafico.currentText()
        modo_spider = "Multielementar" in texto_modo
        norm_condrito = "Condrito" in texto_modo
        em_ppm = ("ppm" in texto_modo) and not norm_condrito

        minerais_validos = KD_DATABASE[magma_atual]
        is_fusao = "Fusão" in processo
        is_afc = "AFC" in processo

        r_afc = self.sl_r.value() / 100.0
        self.lbl_r.setText(f"{r_afc:.2f}")

        ca_c0_global = self.sl_ca_c0.value() / 10.0
        self.lbl_ca_c0.setText(f"{ca_c0_global:.2f}")

        if is_afc:
            ma_50 = (r_afc / (1.0 - r_afc)) * 0.50
            ma_10 = (r_afc / (1.0 - r_afc)) * 0.90
            cor_aviso = "#b35900" if ma_50 > 0.5 else "#333333"
            self.lbl_balanco_afc.setText(
                f"<small><span style='color: {cor_aviso};'><b>Balanço de Massa Assimilada (M<sub>a</sub>/M<sub>m</sub><sup>0</sup>):</b> "
                f"{ma_50:.2f} (em F=0.5) | {ma_10:.2f} (em F=0.1)</span></small>"
            )

        elementos_ativos = []
        for i in range(self.lista_elementos.count()):
            item = self.lista_elementos.item(i)
            if item.checkState() == Qt.Checked:
                elementos_ativos.append(item.text())

        soma_x0 = 0
        soma_p = 0
        bruto_x0 = {}
        bruto_p = {}

        for min_nome in minerais_validos.keys():
            vx = self.sliders_x0[min_nome].value()
            vp = self.sliders_p[min_nome].value()
            self.labels_x0[min_nome].setText(f"{vx}%")
            self.labels_p[min_nome].setText(f"{vp}%")
            bruto_x0[min_nome] = vx
            bruto_p[min_nome] = vp
            soma_x0 += vx
            soma_p += vp

        if soma_x0 == 0:
            self.lbl_soma_x0.setText("<span style='color: red;'><b>Erro: X0 = 0%</b></span>")
            self.limpar_memoria_plot()
            self.ax.clear()
            self.ax.text(0.5, 0.5, "Defina a proporção modal (X0 > 0%) de pelo menos um mineral.", ha='center', va='center', fontsize=11, color='gray')
            self.canvas.draw()
            return
        elif soma_x0 != 100:
            self.lbl_soma_x0.setText(f"<span style='color: orange;'><b>Soma X0: {soma_x0}% (Norm. 100%)</b></span>")
        else:
            self.lbl_soma_x0.setText(f"<span style='color: green;'><b>Soma X0: {soma_x0}%</b></span>")

        if is_fusao:
            if soma_p <= 0:
                self.lbl_soma_p.setText(f"<span style='color: red;'><b>Erro: Soma pi = {soma_p}% (&le; 0%)</b></span>")
                self.limpar_memoria_plot()
                self.ax.clear()
                self.ax.text(
                    0.5, 0.5,
                    "A soma algébrica da reação de fusão (soma de pi) deve ser estritamente positiva (> 0%).\nAjuste pi ou clique em 'Igualar pi = X0'.",
                    ha='center', va='center', fontsize=11, color='gray'
                )
                self.canvas.draw()
                return
            elif soma_p != 100:
                self.lbl_soma_p.setText(f"<span style='color: orange;'><b>Soma pi: {soma_p}% (Norm. 100%)</b></span>")
            else:
                self.lbl_soma_p.setText(f"<span style='color: green;'><b>Soma pi: {soma_p}%</b></span>")

            minerais_impossiveis = [m for m in minerais_validos.keys() if bruto_x0[m] == 0 and bruto_p[m] > 0]
            if minerais_impossiveis:
                nomes_imp = ", ".join(minerais_impossiveis)
                self.lbl_soma_p.setText(f"<span style='color: red;'><b>Erro: {nomes_imp} (X0=0%, pi>0%)</b></span>")
                self.limpar_memoria_plot()
                self.ax.clear()
                self.ax.text(
                    0.5, 0.5,
                    f"Inconsistência Termodinâmica:\nA fase '{nomes_imp}' está sendo consumida na fusão (pi > 0%),\nmas não existe na rocha fonte inicial (X0 = 0%).",
                    ha='center', va='center', fontsize=11, color='#8B0000'
                )
                self.canvas.draw()
                return

        if not elementos_ativos:
            self.limpar_memoria_plot()
            self.ax.clear()
            self.ax.text(0.5, 0.5, "Nenhum elemento selecionado.", ha='center', va='center', fontsize=12, color='gray')
            self.ax.set_title(f"Magma_Trace: {processo}", fontsize=13)
            self.canvas.draw()
            return

        norm_x0 = {k: v / soma_x0 for k, v in bruto_x0.items()}
        norm_p = {k: v / soma_p for k, v in bruto_p.items()} if is_fusao else norm_x0.copy()

        f_max = 1.0
        min_exauridos = []
        if is_fusao:
            for m in minerais_validos.keys():
                if norm_p[m] > norm_x0[m] + 1e-9:
                    f_lim = norm_x0[m] / norm_p[m]
                    if f_lim < f_max - 1e-7:
                        f_max = f_lim
                        min_exauridos = [m]
                    elif abs(f_lim - f_max) <= 1e-7:
                        min_exauridos.append(m)

        bulk_D0 = {el: 0.0 for el in elementos_ativos}
        bulk_P = {el: 0.0 for el in elementos_ativos}
        for el in elementos_ativos:
            for mineral in minerais_validos.keys():
                kd_val = minerais_validos[mineral].get(el, 0.0)
                bulk_D0[el] += norm_x0[mineral] * kd_val
                bulk_P[el] += norm_p[mineral] * kd_val

        F_limite_sup = min(0.995, f_max * 0.998)
        F = np.linspace(0.005, F_limite_sup, 160) if F_limite_sup > 0.005 else np.array([max(1e-4, F_limite_sup)])

        usar_ca_global = "Global" in self.cb_modo_ca.currentText()
        resultados = self.calcular_concentracoes_em_F(
            F, elementos_ativos, bulk_D0, bulk_P, processo, r_afc, ca_c0_global, usar_ca_global,
            em_ppm=em_ppm, norm_condrito=norm_condrito
        )

        self.modo_spider_plot = modo_spider
        self.em_ppm_plot = em_ppm
        self.norm_condrito_plot = norm_condrito
        self.F_plot = F
        self.resultados_plot = resultados
        self.elementos_plot = elementos_ativos
        self.processo_plot = processo
        self.magma_plot = magma_atual
        self.bulk_D0_plot = bulk_D0
        self.bulk_P_plot = bulk_P
        self.moda_x0_plot = norm_x0
        self.moda_p_plot = norm_p
        self.f_max_plot = f_max
        self.min_exauridos_plot = min_exauridos

        self.ax.clear()
        ymin_data = np.inf
        ymax_data = -np.inf

        if not modo_spider:
            self.F_spider_plot = None
            self.res_spider_plot = None
            self.el_ordenados_spider = []

            for el in elementos_ativos:
                cor_idx = ELEMENTOS_MASTER.index(el) % len(CORES_HEX)
                cor = CORES_HEX[cor_idx]
                y_vals = resultados[el]
                ymin_data = min(ymin_data, np.min(y_vals))
                ymax_data = max(ymax_data, np.max(y_vals))

                c0_info = f", $C_0$={C0_DATABASE[el]:.3g} ppm" if em_ppm else ""
                if is_fusao:
                    legenda = f'{el} ($D_0$={bulk_D0[el]:.3f}, $P$={bulk_P[el]:.3f}{c0_info})'
                elif is_afc and not usar_ca_global:
                    raz_real = CA_PPM_DATABASE[el] / max(1e-12, C0_DATABASE[el])
                    legenda = f'{el} ($\\bar{{D}}$={bulk_D0[el]:.3f}, $C_a/C_0$={raz_real:.2g}{c0_info})'
                else:
                    legenda = f'{el} ($\\bar{{D}}$={bulk_D0[el]:.3f}{c0_info})'
                self.ax.plot(F, y_vals, label=legenda, color=cor, linewidth=2.3)

            if is_fusao and min_exauridos and f_max < 0.995:
                nomes_ex = ", ".join(min_exauridos)
                self.ax.axvline(
                    x=f_max, color='#8B0000', linestyle='--', linewidth=1.8,
                    label=f'Exaustão: {nomes_ex} ($F_{{max}}$={f_max:.2f})'
                )

            if is_fusao:
                self.ax.set_xlim(0, 1)
                self.ax.set_xlabel("Fração de Líquido Gerado (F)", fontsize=11)
            else:
                self.ax.set_xlim(1, 0)
                self.ax.set_xlabel("Fração de Líquido Remanescente (F)", fontsize=11)

            self.ax.legend(loc='best', fontsize=8.5, ncol=2 if len(elementos_ativos) > 5 else 1)

        else:
            apenas_etr = all(el in LISTA_ETR for el in elementos_ativos)
            if apenas_etr:
                el_ordenados = [el for el in LISTA_ETR if el in elementos_ativos]
                rotulo_x = "Elementos Terras Raras (Ordem Atômica La -> Lu)"
            else:
                el_ordenados = [el for el in ORDEM_SPIDERGRAM if el in elementos_ativos]
                rotulo_x = "Elementos-Traço (Ordem de Incompatibilidade Manto Primitivo - Sun & McDonough, 1989)"

            if is_fusao:
                if F_limite_sup >= 0.50:
                    F_niveis = np.array([0.01, 0.05, 0.10, 0.25, 0.50])
                elif F_limite_sup > 0.01:
                    F_niveis = np.round(np.linspace(max(0.002, F_limite_sup * 0.05), F_limite_sup, 5), 4)
                    F_niveis = np.unique(F_niveis)
                else:
                    F_niveis = np.array([max(1e-4, F_limite_sup)])
            else:
                F_niveis = np.array([0.90, 0.70, 0.50, 0.30, 0.10])

            res_spider = self.calcular_concentracoes_em_F(
                F_niveis, el_ordenados, bulk_D0, bulk_P, processo, r_afc, ca_c0_global, usar_ca_global,
                em_ppm=em_ppm, norm_condrito=norm_condrito
            )
            res_spider_rel = self.calcular_concentracoes_em_F(
                F_niveis, el_ordenados, bulk_D0, bulk_P, processo, r_afc, ca_c0_global, usar_ca_global,
                em_ppm=False, norm_condrito=False
            )

            self.F_spider_plot = F_niveis
            self.res_spider_plot = res_spider
            self.el_ordenados_spider = el_ordenados

            x_idx = np.arange(len(el_ordenados))
            cores_f = ['#1f77b4', '#2ca02c', '#ff7f0e', '#d62728', '#9467bd']

            for i_f, f_val in enumerate(F_niveis):
                y_spider = np.array([res_spider[el][i_f] for el in el_ordenados])
                ymin_data = min(ymin_data, np.min(y_spider))
                ymax_data = max(ymax_data, np.max(y_spider))

                metricas_extra = []
                if 'La' in el_ordenados and 'Yb' in el_ordenados:
                    if norm_condrito:
                        la_yb_n = res_spider['La'][i_f] / max(1e-30, res_spider['Yb'][i_f])
                        metricas_extra.append(f"(La/Yb)$_N$={la_yb_n:.2f}")
                    elif em_ppm:
                        la_yb_ppm = res_spider['La'][i_f] / max(1e-30, res_spider['Yb'][i_f])
                        metricas_extra.append(f"La/Yb={la_yb_ppm:.2f}")
                    else:
                        la_yb_rel = res_spider_rel['La'][i_f] / max(1e-30, res_spider_rel['Yb'][i_f])
                        metricas_extra.append(f"(La/Yb)$_{{l/0}}$={la_yb_rel:.2f}")

                if all(k in el_ordenados for k in ('Sm', 'Eu', 'Gd')):
                    if norm_condrito:
                        eu_star_n = res_spider['Eu'][i_f] / np.sqrt(max(1e-30, res_spider['Sm'][i_f] * res_spider['Gd'][i_f]))
                        metricas_extra.append(f"Eu/Eu*={eu_star_n:.2f}")
                    elif em_ppm:
                        eu_n = res_spider['Eu'][i_f] / CONDRITO_C1_PPM['Eu']
                        sm_n = res_spider['Sm'][i_f] / CONDRITO_C1_PPM['Sm']
                        gd_n = res_spider['Gd'][i_f] / CONDRITO_C1_PPM['Gd']
                        eu_star_n = eu_n / np.sqrt(max(1e-30, sm_n * gd_n))
                        metricas_extra.append(f"Eu/Eu*={eu_star_n:.2f}")
                    else:
                        eu_star_rel = res_spider_rel['Eu'][i_f] / np.sqrt(max(1e-30, res_spider_rel['Sm'][i_f] * res_spider_rel['Gd'][i_f]))
                        metricas_extra.append(f"(Eu/Eu*)$_{{l/0}}$={eu_star_rel:.2f}")

                sufixo = f" [{', '.join(metricas_extra)}]" if metricas_extra else ""
                self.ax.plot(
                    x_idx, y_spider, marker='o', markersize=5.5, linewidth=2.2,
                    color=cores_f[i_f % len(cores_f)], label=f'F = {f_val:.3g}{sufixo}'
                )

            if norm_condrito:
                y_c0_n = np.array([C0_DATABASE[el] / CONDRITO_C1_PPM[el] for el in el_ordenados])
                ymin_data = min(ymin_data, np.min(y_c0_n))
                ymax_data = max(ymax_data, np.max(y_c0_n))
                self.ax.plot(
                    x_idx, y_c0_n, color='black', linestyle='--', marker='s', markersize=4.5,
                    linewidth=1.3, alpha=0.75, label='Fonte / Magma Inicial ($C_0 / C_{Condrito}$)'
                )
            elif em_ppm:
                y_c0 = np.array([C0_DATABASE[el] for el in el_ordenados])
                ymin_data = min(ymin_data, np.min(y_c0))
                ymax_data = max(ymax_data, np.max(y_c0))
                self.ax.plot(
                    x_idx, y_c0, color='black', linestyle='--', marker='s', markersize=4.5,
                    linewidth=1.3, alpha=0.75, label='Fonte / Magma Inicial ($C_0$ em ppm)'
                )
            else:
                self.ax.axhline(1.0, color='black', linestyle='--', linewidth=1.0, alpha=0.6, label='Fonte / Magma Inicial ($C_0$)')

            self.ax.set_xticks(x_idx)
            self.ax.set_xticklabels(el_ordenados, fontsize=10, fontweight='bold')
            self.ax.set_xlabel(rotulo_x, fontsize=11)
            self.ax.legend(loc='best', fontsize=8.5)

        self.ax.set_yscale('log')

        if self.chk_clamp_y.isChecked() and np.isfinite(ymin_data) and np.isfinite(ymax_data):
            lim_inf_padrao = 1e-5 if (em_ppm or norm_condrito) else 1e-4
            lim_sup_padrao = 1e6 if (em_ppm or norm_condrito) else 1e4
            y_inf = max(lim_inf_padrao, ymin_data * 0.65)
            y_sup = min(lim_sup_padrao, ymax_data * 1.45)
            if y_inf >= y_sup:
                y_inf, y_sup = 1e-2, 1e2
            self.ax.set_ylim(y_inf, y_sup)

        if norm_condrito:
            rotulo_y = "Concentração Normalizada ao Condrito C1 $(C_l / C_{Condrito})$"
        elif em_ppm:
            rotulo_y = "Concentração no Líquido $C_l$ (ppm)"
        else:
            rotulo_y = "Enriquecimento Relativo $(C_l / C_0)$"

        self.ax.set_ylabel(rotulo_y, fontsize=11)
        self.ax.set_title(f"Magma_Trace: {processo}", fontsize=12.5)
        self.ax.grid(True, which="both", linestyle=":", alpha=0.6)

        self.canvas.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    icone_mt, pixmap_mt = construir_icone_mt("magma_trace.ico")
    app.setWindowIcon(icone_mt)

    tela_inicial = TelaBoasVindas(pixmap_icone=pixmap_mt)
    tela_inicial.setWindowIcon(icone_mt)
    if tela_inicial.exec_() == QDialog.Accepted:
        janela = GeochemApp()
        janela.setWindowIcon(icone_mt)
        janela.show()
        sys.exit(app.exec_())