#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Análisis Estático de Marco Plano 2D - Hiperestático
====================================================
Autor: Grupo 4 - Proyecto 1
Descripción: Análisis de marco plano con columna empotrada y viga articulada
             con carga distribuida en columna y carga puntual en viga.
"""

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import sys

# ============================================================================
# 1. GEOMETRÍA Y PROPEDADES MECÁNICAS
# ============================================================================

# Coordenadas de los nodos (metros)
nodos = {
    1: (0.0, 0.0),   # Base empotrada
    2: (0.0, 2.0),   # Unión columna-viga
    3: (0.0, 5.0),   # Extremo superior columna (libre)
    4: (5.0, 2.0),   # Punto de carga puntual en viga
    5: (8.0, 2.0),   # Apoyo derecho (articulado)
}

# Conectividad de elementos [nodo_i, nodo_f]
elementos = {
    1: [1, 2],  # Columna tramo inferior
    2: [2, 3],  # Columna tramo superior
    3: [2, 4],  # Viga tramo izquierdo
    4: [4, 5],  # Viga tramo derecho
}

# Materiales
E_acero = 200e6    # kN/m² (200 GPa)
E_conc = 25e6      # kN/m² (25 GPa - hormigón para viga cuadrada)

# ============================================================================
# 2. CÁLCULO DE PROPIEDADES DE SECCIONES
# ============================================================================

print("=" * 70)
print("PROPIEDADES GEOMÉTRICAS DE LAS SECCIONES")
print("=" * 70)

# --- COLUMNA: Perfil Doble T (I) en Acero ASTM A36 ---
# Dimensiones estándar (en metros)
bf_col = 0.20     # Ancho de ala (mm -> m: 200mm)
tf_col = 0.012    # Espesor de ala (mm -> m: 12mm)
h_col = 0.30      # Alto total (mm -> m: 300mm)
tw_col = 0.008    # Espesor de alma (mm -> m: 8mm)

# Cálculo del área (A = 2*bf*tf + (h - 2*tf)*tw)
A_col = 2 * bf_col * tf_col + (h_col - 2 * tf_col) * tw_col
print(f"\nColumna - Perfil I de Acero ASTM A36:")
print(f"  Ancho de ala (bf) = {bf_col*1000:.0f} mm")
print(f"  Espesor de ala (tf) = {tf_col*1000:.0f} mm")
print(f"  Alto total (h) = {h_col*1000:.0f} mm")
print(f"  Espesor de alma (tw) = {tw_col*1000:.0f} mm")
print(f"  Área (A_col) = {A_col:.6f} m² = {A_col*1e4:.2f} cm²")

# Momento de inercia respecto al eje fuerte (Ix)
# Ix = (bf*h³/12) - ((bf-tw)*(h-2*tf)³/12)
I_col = (bf_col * h_col**3 / 12) - ((bf_col - tw_col) * (h_col - 2*tf_col)**3 / 12)
print(f"  Momento de Inercia (I_col) = {I_col:.8f} m⁴ = {I_col*1e8:.2f} cm⁴")

# --- VIGA: Sección cuadrada de 40 cm x 40 cm (Hormigón) ---
b_viga = 0.40     # Ancho (m)
h_viga = 0.40     # Alto (m)

A_viga = b_viga * h_viga
I_viga = (b_viga * h_viga**3) / 12

print(f"\nViga - Sección cuadrada de Hormigón:")
print(f"  Dimensión = {b_viga*100:.0f} cm x {h_viga*100:.0f} cm")
print(f"  Área (A_viga) = {A_viga:.4f} m² = {A_viga*1e4:.0f} cm²")
print(f"  Momento de Inercia (I_viga) = {I_viga:.8f} m⁴ = {I_viga*1e8:.2f} cm⁴")

# ============================================================================
# 3. ANÁLISIS CON OPENSEESPY
# ============================================================================

print("\n" + "=" * 70)
print("ANÁLISIS ESTRUCTURAL CON OPENSEESPY")
print("=" * 70)

# Limpiar el modelo
ops.wipe()

# Crear modelo (2D, 3 grados de libertad por nodo: UX, UY, RZ)
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- Definir material elástico ---
matTag_col = 1
matTag_viga = 2
ops.uniaxialMaterial('Elastic', matTag_col, E_acero)
ops.uniaxialMaterial('Elastic', matTag_viga, E_conc)

# --- Definir nodos ---
for nodo_id, (x, y) in nodos.items():
    ops.node(nodo_id, x, y)

# --- Definir secciones geométricas ---
# Columna: sección I
ops.section('Elastic', 1, E_acero, A_col, I_col, 1.0, 1.0, 1.0)

# Viga: sección cuadrada
ops.section('Elastic', 2, E_conc, A_viga, I_viga, 1.0, 1.0, 1.0)

# --- Definir transformación geométrica ---
ops.geomTransf('Linear', 1)   # Para columna (vertical)
ops.geomTransf('Linear', 2)   # Para viga (horizontal)

# --- Definir elementos ---
# Columna
ops.element('elasticBeamColumn', 1, 1, 2, A_col, E_acero, I_col, 1)
ops.element('elasticBeamColumn', 2, 2, 3, A_col, E_acero, I_col, 1)

# Viga
ops.element('elasticBeamColumn', 3, 2, 4, A_viga, E_conc, I_viga, 2)
ops.element('elasticBeamColumn', 4, 4, 5, A_viga, E_conc, I_viga, 2)

# --- Condiciones de borde ---
# Nodo 1: Empotrado (restringe todos los DOF)
ops.fix(1, 1, 1, 1)

# Nodo 5: Articulado (restringe UX y UY, permite rotación)
ops.fix(5, 1, 1, 0)

# --- Cargas ---
# Tiempo de carga
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

# Carga distribuida en columna (q = 17 kN/m en dirección +X)
# Se aplica como fuerza equivalente en los nodos
q = 17.0  # kN/m

# Columna elemental 1 (Nodo 1 a 2, altura = 2 m)
L1 = 2.0
Fx1 = q * L1 / 2  # Fuerza en Nodo 1
Fx2 = q * L1 / 2  # Fuerza en Nodo 2
ops.load(1, Fx1, 0.0, 0.0)
ops.load(2, Fx2, 0.0, 0.0)

# Columna elemental 2 (Nodo 2 a 3, altura = 3 m)
L2 = 3.0
Fx2b = q * L2 / 2  # Fuerza adicional en Nodo 2
Fx3 = q * L2 / 2    # Fuerza en Nodo 3
ops.load(2, Fx2b, 0.0, 0.0)  # Se suma a la carga anterior
ops.load(3, Fx3, 0.0, 0.0)

# Carga puntual en Nodo 4 (F = 20 kN en dirección -Y)
ops.load(4, 0.0, -20.0, 0.0)

# --- Análisis ---
ops.constraints('Plain')
ops.numberer('RCM')
ops.system('BandGeneral')
ops.test('NormDispIncr', 1e-10, 10)
ops.algorithm('Linear')
ops.integrator('LoadControl', 1)
ops.algorithm('Linear')
ops.analysis('Static')
ok = ops.analyze(1)

if ok == 0:
    print("\n¡Análisis completado exitosamente!")
else:
    print("\nError en el análisis")
    sys.exit(1)

# ============================================================================
# 4. RESULTADOS - REACCIONES DE APOYO
# ============================================================================

print("\n" + "=" * 70)
print("REACCIONES DE APOYO")
print("=" * 70)

# Reacciones en Nodo 1 (Empotramiento)
rxn1 = ops.nodeReaction(1)
print(f"\nNodo 1 (Empotrado):")
print(f"  Fx = {rxn1[0]:.4f} kN")
print(f"  Fy = {rxn1[1]:.4f} kN")
print(f"  Mz = {rxn1[2]:.4f} kN·m")

# Reacciones en Nodo 5 (Articulado)
rxn5 = ops.nodeReaction(5)
print(f"\nNodo 5 (Articulado):")
print(f"  Fx = {rxn5[0]:.4f} kN")
print(f"  Fy = {rxn5[1]:.4f} kN")
print(f"  Mz = {rxn5[2]:.4f} kN·m")

# Verificación de equilibrio
Fx_total_carga = q * 5  # Carga distribuida total
Fy_total_carga = -20.0  # Carga puntual

Fx_reacciones = rxn1[0] + rxn5[0]
Fy_reacciones = rxn1[1] + rxn5[1]
M_reacciones = rxn1[2] + rxn5[2]

print(f"\n--- Verificación de Equilibrio ---")
print(f"ΣFx = {Fx_reacciones:.4f} kN (carga aplicada: {Fx_total_carga:.4f} kN)")
print(f"ΣFy = {Fy_reacciones:.4f} kN (carga aplicada: {Fy_total_carga:.4f} kN)")

# ============================================================================
# 5. DESPLAZAMIENTOS DE NODOS
# ============================================================================

print("\n" + "=" * 70)
print("DESPLAZAMIENTOS DE NODOS")
print("=" * 70)

for nodo_id in nodos:
    disp = ops.nodeDisp(nodo_id)
    print(f"Nodo {nodo_id}: UX = {disp[0]*1000:.4f} mm, "
          f"UY = {disp[1]*1000:.4f} mm, "
          f"RZ = {np.degrees(disp[2]):.6f}°")

# ============================================================================
# 6. FUERZAS INTERNAS POR ELEMENTO
# ============================================================================

print("\n" + "=" * 70)
print("FUERZAS INTERNAS POR ELEMENTO")
print("=" * 70)

for elem_id in range(1, 5):
    fuerzas = ops.eleForce(elem_id)
    # Formato: [N_i, V_i, M_i, N_j, V_j, M_j]
    print(f"\nElemento {elem_id} (Nodos {elementos[elem_id][0]}-{elementos[elem_id][1]}):")
    print(f"  Nodo i: N = {fuerzas[0]:.4f} kN, V = {fuerzas[1]:.4f} kN, M = {fuerzas[2]:.4f} kN·m")
    print(f"  Nodo j: N = {fuerzas[3]:.4f} kN, V = {fuerzas[4]:.4f} kN, M = {fuerzas[5]:.4f} kN·m")

# ============================================================================
# 7. DIAGRAMAS DE FUERZAS INTERNAS
# ============================================================================

print("\n" + "=" * 70)
print("DIAGRAMAS DE FUERZAS INTERNAS")
print("=" * 70)

fig, axes = plt.subplots(3, 1, figsize=(12, 14))
fig.suptitle('Diagramas de Fuerzas Internas - Marco Plano 2D', fontsize=14, fontweight='bold')

# Colores para los elementos
colores = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']
nombres_elem = ['Columna 1-2', 'Columna 2-3', 'Viga 2-4', 'Viga 4-5']

# Para cada elemento, obtener fuerzas internas y graficar
for elem_id in range(1, 5):
    fuerzas = ops.eleForce(elem_id)
    nodo_i, nodo_f = elementos[elem_id]
    xi, yi = nodos[nodo_i]
    xf, yf = nodos[nodo_f]
    
    # Longitud del elemento
    L_elem = np.sqrt((xf - xi)**2 + (yf - yi)**2)
    
    # Puntos a lo largo del elemento (para graficar)
    n_pts = 20
    t_vals = np.linspace(0, 1, n_pts)
    
    # Coordenadas x para graficar (proyectadas)
    x_vals = xi + t_vals * (xf - xi)
    
    # --- Diagrama de Cortante (V) ---
    # Distribución lineal de cortante
    V_i = fuerzas[1]
    V_j = fuerzas[4]
    V_vals = V_i + (V_j - V_i) * t_vals
    
    # Para la columna con carga distribuida, ajustar el cortante
    if elem_id in [1, 2]:  # Elementos de columna
        q_elem = q  # 17 kN/m
        # El cortante varía linealmente con la carga distribuida
        V_vals = V_i + q_elem * (xf - xi) * t_vals
    
    axes[0].plot(x_vals, V_vals, color=colores[elem_id-1], linewidth=2, label=nombres_elem[elem_id-1])
    axes[0].fill_between(x_vals, V_vals, alpha=0.1, color=colores[elem_id-1])
    
    # --- Diagrama de Momento Flector (M) ---
    # Distribución cuadrática del momento para carga distribuida
    M_i = fuerzas[2]
    M_j = fuerzas[5]
    
    if elem_id in [1, 2]:  # Elementos de columna con carga distribuida
        # Momento parabólico: M(t) = M_i + V_i*t*L + q*t²*L²/2
        M_vals = M_i + V_i * t_vals * L_elem + q_elem * (t_vals * L_elem)**2 / 2
    else:  # Viga sin carga distribuida (solo puntual)
        # Momento lineal entre nodos (sin carga distribuida en viga)
        M_vals = M_i + (M_j - M_i) * t_vals
    
    axes[2].plot(x_vals, M_vals, color=colores[elem_id-1], linewidth=2, label=nombres_elem[elem_id-1])
    axes[2].fill_between(x_vals, M_vals, alpha=0.1, color=colores[elem_id-1])

# --- Diagrama de Axial (N) ---
# Para la columna con carga horizontal, hay fuerza axial por deformación
# Simplificación: mostrar fuerza axial constante por elemento
for elem_id in range(1, 5):
    fuerzas = ops.eleForce(elem_id)
    nodo_i, nodo_f = elementos[elem_id]
    xi, yi = nodos[nodo_i]
    xf, yf = nodos[nodo_f]
    
    x_vals = np.linspace(xi, xf, 20)
    N_i = fuerzas[0]
    N_j = fuerzas[3]
    N_vals = np.full_like(x_vals, N_i)  # Axial constante (simplificación)
    
    axes[1].plot(x_vals, N_vals, color=colores[elem_id-1], linewidth=2, label=nombres_elem[elem_id-1])
    axes[1].fill_between(x_vals, N_vals, alpha=0.1, color=colores[elem_id-1])

# Configurar ejes
axes[0].set_title('Diagrama de Cortante (V) [kN]', fontsize=12)
axes[0].set_xlabel('Posición X [m]')
axes[0].set_ylabel('Cortante [kN]')
axes[0].grid(True, alpha=0.3)
axes[0].axhline(y=0, color='k', linewidth=0.5)
axes[0].legend()

axes[1].set_title('Diagrama de Fuerza Axial (N) [kN]', fontsize=12)
axes[1].set_xlabel('Posición X [m]')
axes[1].set_ylabel('Axial [kN]')
axes[1].grid(True, alpha=0.3)
axes[1].axhline(y=0, color='k', linewidth=0.5)
axes[1].legend()

axes[2].set_title('Diagrama de Momento Flector (M) [kN·m]', fontsize=12)
axes[2].set_xlabel('Posición X [m]')
axes[2].set_ylabel('Momento [kN·m]')
axes[2].grid(True, alpha=0.3)
axes[2].axhline(y=0, color='k', linewidth=0.5)
axes[2].legend()

plt.tight_layout()
plt.savefig('diagramas_fuerzas_internas.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nDiagramas guardados como 'diagramas_fuerzas_internas.png'")

# ============================================================================
# 8. RESUMEN DE RESULTADOS
# ============================================================================

print("\n" + "=" * 70)
print("RESUMEN GENERAL DE RESULTADOS")
print("=" * 70)

print(f"""
┌─────────────────────────────────────────────────────────────┐
│                   PROPIEDADES DE SECCIONES                  │
├─────────────────────────────────────────────────────────────┤
│ Columna (Acero A36):                                        │
│   A = {A_col:.6f} m² ({A_col*1e4:.2f} cm²)                          │
│   I = {I_col:.8f} m⁴ ({I_col*1e8:.2f} cm⁴)                   │
│                                                             │
│ Viga (Hormigón 40x40):                                      │
│   A = {A_viga:.4f} m² ({A_viga*1e4:.0f} cm²)                            │
│   I = {I_viga:.8f} m⁴ ({I_viga*1e8:.2f} cm⁴)                    │
├─────────────────────────────────────────────────────────────┤
│                    REACCIONES DE APOYO                       │
├─────────────────────────────────────────────────────────────┤
│ Nodo 1 (Empotrado):                                         │
│   Fx = {rxn1[0]:.4f} kN                                            │
│   Fy = {rxn1[1]:.4f} kN                                            │
│   Mz = {rxn1[2]:.4f} kN·m                                         │
│                                                             │
│ Nodo 5 (Articulado):                                        │
│   Fx = {rxn5[0]:.4f} kN                                            │
│   Fy = {rxn5[1]:.4f} kN                                            │
├─────────────────────────────────────────────────────────────┤
│                  CARGAS APLICADAS                            │
├─────────────────────────────────────────────────────────────┤
│ Columna: q = 17 kN/m (dirección +X, altura total 5m)        │
│ Viga: F = 20 kN (dirección -Y, en X = 5m)                   │
└─────────────────────────────────────────────────────────────┘
""")

# ============================================================================
# 9. EXPORTAR RESULTADOS A ARCHIVO
# ============================================================================

with open('resultados_analisis.txt', 'w', encoding='utf-8') as f:
    f.write("RESULTADOS DEL ANÁLISIS ESTÁTICO - MARCO PLANO 2D\n")
    f.write("=" * 60 + "\n\n")
    
    f.write("PROPIEDADES DE SECCIONES:\n")
    f.write(f"Columna: A = {A_col:.6f} m², I = {I_col:.8f} m⁴\n")
    f.write(f"Viga: A = {A_viga:.4f} m², I = {I_viga:.8f} m⁴\n\n")
    
    f.write("REACCIONES DE APOYO:\n")
    f.write(f"Nodo 1 (Empotrado): Fx = {rxn1[0]:.4f} kN, Fy = {rxn1[1]:.4f} kN, Mz = {rxn1[2]:.4f} kN·m\n")
    f.write(f"Nodo 5 (Articulado): Fx = {rxn5[0]:.4f} kN, Fy = {rxn5[1]:.4f} kN\n\n")
    
    f.write("FUERZAS INTERNAS:\n")
    for elem_id in range(1, 5):
        fuerzas = ops.eleForce(elem_id)
        f.write(f"Elemento {elem_id}: N_i={fuerzas[0]:.4f}, V_i={fuerzas[1]:.4f}, M_i={fuerzas[2]:.4f}\n")
        f.write(f"              N_j={fuerzas[3]:.4f}, V_j={fuerzas[4]:.4f}, M_j={fuerzas[5]:.4f}\n")

print("Resultados exportados a 'resultados_analisis.txt'")

# Limpiar modelo
ops.wipe()
print("\n¡Análisis completado!")
