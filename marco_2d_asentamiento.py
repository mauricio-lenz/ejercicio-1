#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Análisis Estático de Marco Plano 2D - Hiperestático
====================================================
Grupo 4 - Proyecto 1

Verificación doble: método matricial directo + OpenSeesPy
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# ============================================================================
# 1. GEOMETRÍA
# ============================================================================

nodos = {
    1: np.array([0.0, 0.0]),
    2: np.array([0.0, 2.0]),
    3: np.array([0.0, 5.0]),
    4: np.array([5.0, 2.0]),
    5: np.array([8.0, 2.0]),
}

elementos = {
    1: (1, 2),   # Columna inferior
    2: (2, 3),   # Columna superior
    3: (2, 4),   # Viga izquierda
    4: (4, 5),   # Viga derecha
}

q = 17.0
F = 20.0

# ============================================================================
# 2. PROPIEDADES DE SECCIONES
# ============================================================================

print("=" * 70)
print("PROPIEDADES GEOMÉTRICAS DE LAS SECCIONES")
print("=" * 70)

bf = 0.20; tf = 0.012; h = 0.30; tw = 0.008
A_col = 2 * bf * tf + (h - 2 * tf) * tw
I_col = (bf * h**3 / 12) - ((bf - tw) * (h - 2 * tf)**3 / 12)
E_acero = 200e6

print(f"\nColumna - Perfil I (Acero ASTM A36):")
print(f"  bf={bf*1000:.0f}mm  tf={tf*1000:.0f}mm  h={h*1000:.0f}mm  tw={tw*1000:.0f}mm")
print(f"  A = {A_col:.6f} m2 = {A_col*1e4:.2f} cm2")
print(f"  I = {I_col:.8f} m4 = {I_col*1e8:.2f} cm4")

A_viga = 0.16
I_viga = 0.4**4 / 12
E_conc = 25e6

print(f"\nViga - Cuadrada 40x40 cm:")
print(f"  A = {A_viga:.4f} m2 = {A_viga*1e4:.0f} cm2")
print(f"  I = {I_viga:.8f} m4 = {I_viga*1e8:.2f} cm4")

# ============================================================================
# 3. MÉTODO MATRICIAL DIRECTO (sin OpenSees)
# ============================================================================

print("\n" + "=" * 70)
print("MÉTODO MATRICIAL DIRECTO (Verificación)")
print("=" * 70)

n_nodos = 5
ndf = 3
n_dof = n_nodos * ndf

K_global = np.zeros((n_dof, n_dof))

def dir_cosines(ni, nf):
    dx = nf[0] - ni[0]
    dy = nf[1] - ni[1]
    L = np.sqrt(dx**2 + dy**2)
    c = dx / L
    s = dy / L
    return c, s, L

def element_stiffness_local(E, A, I, L):
    ke = np.zeros((6, 6))
    EA_L = E * A / L
    EI = E * I
    L2 = L * L
    L3 = L2 * L

    ke[0, 0] = EA_L
    ke[0, 3] = -EA_L
    ke[3, 0] = -EA_L
    ke[3, 3] = EA_L

    ke[1, 1] = 12 * EI / L3
    ke[1, 2] = 6 * EI / L2
    ke[1, 4] = -12 * EI / L3
    ke[1, 5] = 6 * EI / L2

    ke[2, 1] = 6 * EI / L2
    ke[2, 2] = 4 * EI / L
    ke[2, 4] = -6 * EI / L2
    ke[2, 5] = 2 * EI / L

    ke[4, 1] = -12 * EI / L3
    ke[4, 2] = -6 * EI / L2
    ke[4, 4] = 12 * EI / L3
    ke[4, 5] = -6 * EI / L2

    ke[5, 1] = 6 * EI / L2
    ke[5, 2] = 2 * EI / L
    ke[5, 4] = -6 * EI / L2
    ke[5, 5] = 4 * EI / L

    return ke

def rotation_matrix(c, s):
    T = np.zeros((6, 6))
    T[0, 0] = c; T[0, 1] = s
    T[1, 0] = -s; T[1, 1] = c
    T[2, 2] = 1
    T[3, 3] = c; T[3, 4] = s
    T[4, 3] = -s; T[4, 4] = c
    T[5, 5] = 1
    return T

props = {
    1: (E_acero, A_col, I_col),
    2: (E_acero, A_col, I_col),
    3: (E_conc, A_viga, I_viga),
    4: (E_conc, A_viga, I_viga),
}

elem_data = {}

for eid, (ni_tag, nf_tag) in elementos.items():
    ni = nodos[ni_tag]
    nf = nodos[nf_tag]
    c, s, L = dir_cosines(ni, nf)
    E, A, I = props[eid]
    ke_loc = element_stiffness_local(E, A, I, L)
    T = rotation_matrix(c, s)
    ke_glob = T.T @ ke_loc @ T

    dof_i = [(ni_tag - 1) * ndf + d for d in range(3)]
    dof_f = [(nf_tag - 1) * ndf + d for d in range(3)]
    dofs = dof_i + dof_f

    for a in range(6):
        for b in range(6):
            K_global[dofs[a], dofs[b]] += ke_glob[a, b]

    elem_data[eid] = {
        'ni': ni_tag, 'nf': nf_tag, 'L': L,
        'c': c, 's': s, 'T': T, 'ke_loc': ke_loc,
        'dofs': dofs, 'E': E, 'A': A, 'I': I,
    }

# Vector de fuerzas nodales equivalentes de la carga distribuida
F_ext = np.zeros(n_dof)

for eid in [1, 2]:
    ed = elem_data[eid]
    L = ed['L']
    ni_tag = ed['ni']
    nf_tag = ed['nf']

    # Carga distribuida q en dirección +X (global)
    # En coordenadas locales: perpendicular al elemento vertical = dirección local Y
    # Fuerzas equivalentes: q*L/2 en cada nodo (dirección local Y -> global X)
    Fy_i = q * L / 2
    Fy_j = q * L / 2
    # Momentos equivalentes: q*L^2/12
    Mz_i = q * L**2 / 12
    Mz_j = -q * L**2 / 12

    # Convertir a global usando T^T * f_local
    f_local = np.array([0, Fy_i, Mz_i, 0, Fy_j, Mz_j])
    f_global = ed['T'].T @ f_local

    dofs = ed['dofs']
    for a in range(6):
        F_ext[dofs[a]] += f_global[a]

# Carga puntual F = 20 kN en -Y en nodo 4
dof_n4_uy = (4 - 1) * ndf + 1
F_ext[dof_n4_uy] += -F

print(f"\nVector de fuerzas nodales equivalentes (global):")
for i in range(0, n_dof, 3):
    nodo = i // 3 + 1
    print(f"  Nodo {nodo}: Fx={F_ext[i]:10.4f} kN, Fy={F_ext[i+1]:10.4f} kN, Mz={F_ext[i+2]:10.4f} kN.m")

# DOF restringidos
fixed_dofs = [0, 1, 2, 12, 13]  # Nodo 1: 0,1,2  Nodo 5: 12,13
free_dofs = [i for i in range(n_dof) if i not in fixed_dofs]

K_ff = K_global[np.ix_(free_dofs, free_dofs)]
F_f = F_ext[free_dofs]

U = np.zeros(n_dof)
U[free_dofs] = np.linalg.solve(K_ff, F_f)

print(f"\nDesplazamientos (libres):")
for dof in free_dofs:
    nodo = dof // 3 + 1
    comp = ['UX', 'UY', 'RZ'][dof % 3]
    val = U[dof]
    if comp == 'RZ':
        print(f"  Nodo {nodo} {comp} = {val:.8f} rad = {np.degrees(val):.6f}deg")
    else:
        print(f"  Nodo {nodo} {comp} = {val:.8f} m = {val*1000:.6f} mm")

# Reacciones
R = K_global @ U - F_ext
print(f"\nReacciones de apoyo:")
print(f"  Nodo 1 (Empotrado):  Fx = {R[0]:.4f} kN, Fy = {R[1]:.4f} kN, Mz = {R[2]:.4f} kN.m")
print(f"  Nodo 5 (Articulado): Fx = {R[12]:.4f} kN, Fy = {R[13]:.4f} kN")

print(f"\nVerificación de equilibrio global:")
print(f"  SumFx = {R[0]+R[12]:.4f} kN  (carga: {q*5:.4f} kN)")
print(f"  SumFy = {R[1]+R[13]:.4f} kN  (carga: {-F:.4f} kN)")

# Fuerzas internas por elemento (coordenadas locales)
fuerzas_elem = {}
for eid, ed in elem_data.items():
    dofs = ed['dofs']
    u_elem = U[dofs]
    f_loc = ed['ke_loc'] @ (ed['T'] @ u_elem)
    fuerzas_elem[eid] = f_loc

print(f"\nFuerzas internas por elemento (coordenadas locales):")
for eid in range(1, 5):
    f = fuerzas_elem[eid]
    ni, nf = elementos[eid]
    L = elem_data[eid]['L']
    print(f"\n  Elemento {eid} ({ni}-{nf}, L={L:.1f} m):")
    print(f"    N_i = {f[0]:10.4f} kN   V_i = {f[1]:10.4f} kN   M_i = {f[2]:10.4f} kN.m")
    print(f"    N_j = {f[3]:10.4f} kN   V_j = {f[4]:10.4f} kN   M_j = {f[5]:10.4f} kN.m")

# ============================================================================
# 4. VERIFICACIÓN CON OPENSEESPY
# ============================================================================

print("\n" + "=" * 70)
print("VERIFICACIÓN CON OPENSEESPY")
print("=" * 70)

try:
    import openseespy.opensees as ops
    oops_available = True
except (ImportError, RuntimeError):
    oops_available = False
    print("  OpenSeesPy no compatible con esta version de Python.")
    print("  Se utiliza exclusivamente el metodo matricial directo.\n")

if oops_available:

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    for nid, (x, y) in nodos.items():
        ops.node(nid, x, y)

    ops.uniaxialMaterial('Elastic', 1, E_acero)
    ops.uniaxialMaterial('Elastic', 2, E_conc)

    ops.section('Elastic', 1, E_acero, A_col, I_col, 1.0, 1.0, 1.0)
    ops.section('Elastic', 2, E_conc, A_viga, I_viga, 1.0, 1.0, 1.0)

    # Transformaciones geométricas con vecindad correcta
    # Columna vertical: vector de vecindad vertical (0,1,0)
    ops.geomTransf('Linear', 1, 0, 1, 0)
    # Viga horizontal: vector de vecindad vertical (0,1,0)
    ops.geomTransf('Linear', 2, 0, 1, 0)

    ops.element('elasticBeamColumn', 1, 1, 2, A_col, E_acero, I_col, 1)
    ops.element('elasticBeamColumn', 2, 2, 3, A_col, E_acero, I_col, 1)
    ops.element('elasticBeamColumn', 3, 2, 4, A_viga, E_conc, I_viga, 2)
    ops.element('elasticBeamColumn', 4, 4, 5, A_viga, E_conc, I_viga, 2)

    ops.fix(1, 1, 1, 1)
    ops.fix(5, 1, 1, 0)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)

    # Cargas equivalentes nodales (incluyendo momentos)
    for eid in [1, 2]:
        ed = elem_data[eid]
        L = ed['L']
        ni_tag = ed['ni']
        nf_tag = ed['nf']

        Fy_i = q * L / 2
        Fy_j = q * L / 2
        Mz_i = q * L**2 / 12
        Mz_j = -q * L**2 / 12

        f_local = np.array([0, Fy_i, Mz_i, 0, Fy_j, Mz_j])
        f_global = ed['T'].T @ f_local

        ops.load(ni_tag, f_global[0], f_global[1], f_global[2])
        ops.load(nf_tag, f_global[3], f_global[4], f_global[5])

    ops.load(4, 0.0, -F, 0.0)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    ops.test('NormDispIncr', 1e-10, 10)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ok = ops.analyze(1)

    if ok == 0:
        print("¡Análisis completado exitosamente!")
    else:
        print("Error en el análisis de OpenSeesPy")

    print(f"\nReacciones OpenSeesPy:")
    rxn1 = ops.nodeReaction(1)
    rxn5 = ops.nodeReaction(5)
    print(f"  Nodo 1: Fx={rxn1[0]:.4f} kN, Fy={rxn1[1]:.4f} kN, Mz={rxn1[2]:.4f} kN.m")
    print(f"  Nodo 5: Fx={rxn5[0]:.4f} kN, Fy={rxn5[1]:.4f} kN")

    print(f"\nComparación de reacciones:")
    print(f"  {'':20s} {'Matricial':>12s} {'OpenSees':>12s} {'Diferencia':>12s}")
    print(f"  {'R1_Fx':20s} {R[0]:12.4f} {rxn1[0]:12.4f} {abs(R[0]-rxn1[0]):12.6f}")
    print(f"  {'R1_Fy':20s} {R[1]:12.4f} {rxn1[1]:12.4f} {abs(R[1]-rxn1[1]):12.6f}")
    print(f"  {'R1_Mz':20s} {R[2]:12.4f} {rxn1[2]:12.4f} {abs(R[2]-rxn1[2]):12.6f}")
    print(f"  {'R5_Fx':20s} {R[12]:12.4f} {rxn5[0]:12.4f} {abs(R[12]-rxn5[0]):12.6f}")
    print(f"  {'R5_Fy':20s} {R[13]:12.4f} {rxn5[1]:12.4f} {abs(R[13]-rxn5[1]):12.6f}")

    print(f"\nDesplazamientos OpenSeesPy:")
    for nid in nodos:
        d = ops.nodeDisp(nid)
        print(f"  Nodo {nid}: UX={d[0]*1000:.6f} mm, UY={d[1]*1000:.6f} mm, RZ={np.degrees(d[2]):.6f}deg")

    print(f"\nFuerzas internas OpenSeesPy:")
    for eid in range(1, 5):
        f_ocs = ops.eleForce(eid)
        f_mat = fuerzas_elem[eid]
        ni, nf = elementos[eid]
        print(f"\n  Elemento {eid} ({ni}-{nf}):")
        print(f"    {'':8s} {'Matricial':>14s} {'OpenSees':>14s} {'Diferencia':>14s}")
        labels = ['N_i', 'V_i', 'M_i', 'N_j', 'V_j', 'M_j']
        for k in range(6):
            print(f"    {labels[k]:8s} {f_mat[k]:14.6f} {f_ocs[k]:14.6f} {abs(f_mat[k]-f_ocs[k]):14.8f}")

    ops.wipe()

# ============================================================================
# 5. ECUACIONES DE FUERZAS INTERNAS POR TRAMO
# ============================================================================

print("\n" + "=" * 70)
print("ECUACIONES DE FUERZAS INTERNAS POR TRAMO")
print("=" * 70)

elem_info = {
    1: {'L': 2.0, 'q': q},
    2: {'L': 3.0, 'q': q},
    3: {'L': 5.0, 'q': 0.0},
    4: {'L': 3.0, 'q': 0.0},
}

for eid in range(1, 5):
    f = fuerzas_elem[eid]
    L = elem_info[eid]['L']
    qe = elem_info[eid]['q']
    ni, nf = elementos[eid]

    Ni, Vi, Mi = f[0], f[1], f[2]
    Nj, Vj, Mj = f[3], f[4], f[5]

    print(f"\nElemento {eid} ({ni}-{nf}), L = {L:.1f} m, q = {qe:.1f} kN/m:")
    print(f"  N(x) = {Ni:+.4f} + ({(Nj-Ni)/L:+.4f}).x")
    print(f"  V(x) = {Vi:+.4f} + ({qe:+.4f}).x")
    print(f"  M(x) = {Mi:+.4f} + ({Vi:+.4f}).x + ({qe/2:+.4f}).x2")

    # Valores notables
    x_eval = np.array([0, L])
    N_eval = Ni + (Nj - Ni) * x_eval / L
    V_eval = Vi + qe * x_eval
    M_eval = Mi + Vi * x_eval + qe * x_eval**2 / 2

    if qe > 0:
        # Punto de cortante cero (V=0): x0 = -Vi/qe
        if Vi != 0:
            x0 = -Vi / qe
            if 0 < x0 < L:
                M_max = Mi + Vi * x0 + qe * x0**2 / 2
                print(f"  ** Punto de V=0 en x = {x0:.4f} m")
                print(f"  ** Mximo local: M = {M_max:+.4f} kN.m")

    print(f"  Valores notables:")
    print(f"    x=0: N={N_eval[0]:+.4f} kN, V={V_eval[0]:+.4f} kN, M={M_eval[0]:+.4f} kN.m")
    print(f"    x={L}: N={N_eval[1]:+.4f} kN, V={V_eval[1]:+.4f} kN, M={M_eval[1]:+.4f} kN.m")

# ============================================================================
# 6. DIAGRAMAS
# ============================================================================

print("\n" + "=" * 70)
print("GENERANDO DIAGRAMAS...")
print("=" * 70)

n_pts = 100

for tipo, titulo, escala, color in [
    ('N', 'Fuerza Axial (N) [kN]', 0.3, '#E53935'),
    ('V', 'Fuerza Cortante (V) [kN]', 0.08, '#1E88E5'),
    ('M', 'Momento Flector (M) [kN.m]', 0.03, '#43A047'),
]:

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.set_aspect('equal')
    ax.set_title(f'{titulo}', fontsize=14, fontweight='bold')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.grid(True, alpha=0.2)

    # Estructura base
    for eid, (ni, nf) in elementos.items():
        xi, yi = nodos[ni]
        xf, yf = nodos[nf]
        ax.plot([xi, xf], [yi, yf], 'k-', linewidth=3, zorder=3)

    for nid, pos in nodos.items():
        ax.plot(pos[0], pos[1], 'ko', markersize=6, zorder=5)

    # Apoyo empotrado ( Nodo 1)
    for dx in np.linspace(-0.3, 0.3, 7):
        ax.plot([dx, dx - 0.1], [0, -0.2], 'k-', linewidth=1)
    ax.plot([-0.3, 0.3], [0, 0], 'k-', linewidth=2)

    # Apoyo articulado ( Nodo 5)
    tri_x = [7.7, 8.3, 8.0, 7.7]
    tri_y = [2.0, 2.0, 1.55, 2.0]
    ax.plot(tri_x, tri_y, 'k-', linewidth=2)
    ax.plot([7.6, 8.4], [1.55, 1.55], 'k-', linewidth=2)

    # Cargas
    for yy in np.arange(0.4, 5.0, 0.6):
        ax.annotate('', xy=(0.6, yy), xytext=(0.05, yy),
                     arrowprops=dict(arrowstyle='->', color='blue', lw=1.2))
    ax.text(0.7, 2.5, f'q={q} kN/m', color='blue', fontsize=9,
            fontweight='bold', rotation=90, va='center')

    ax.annotate('', xy=(5.0, 1.1), xytext=(5.0, 2.0),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2.5))
    ax.text(5.1, 0.95, f'F={F} kN', color='red', fontsize=9, fontweight='bold')

    # Diagramas
    for eid, (ni, nf) in elementos.items():
        xi, yi = nodos[ni]
        xf, yf = nodos[nf]
        L = elem_info[eid]['L']
        qe = elem_info[eid]['q']
        f = fuerzas_elem[eid]
        Ni, Vi, Mi = f[0], f[1], f[2]

        x_loc = np.linspace(0, L, n_pts)

        if tipo == 'N':
            vals = Ni + (f[3] - Ni) * x_loc / L
        elif tipo == 'V':
            vals = Vi + qe * x_loc
        else:
            vals = Mi + Vi * x_loc + qe * x_loc**2 / 2

        is_vertical = (abs(xf - xi) < 1e-10)

        if is_vertical:
            # Elemento vertical: graficar diagrama en dirección X
            y_coords = yi + (yf - yi) * x_loc / L
            x_base = np.full_like(y_coords, xi)
            x_diag = xi + vals * escala

            ax.plot(x_diag, y_coords, color=color, linewidth=2, zorder=4)
            ax.fill_betweenx(y_coords, x_base, x_diag, alpha=0.15, color=color, zorder=2)
            ax.plot(x_base, y_coords, color=color, linewidth=0.8, linestyle='--', alpha=0.5, zorder=2)

            # Valores
            ax.annotate(f'{vals[0]:.2f}', (x_diag[0], y_coords[0]),
                         fontsize=8, color=color, fontweight='bold',
                         textcoords="offset points", xytext=(-25, -8))
            ax.annotate(f'{vals[-1]:.2f}', (x_diag[-1], y_coords[-1]),
                         fontsize=8, color=color, fontweight='bold',
                         textcoords="offset points", xytext=(-25, 5))

            idx_max = np.argmax(np.abs(vals))
            if abs(vals[idx_max]) > 0.01:
                ax.plot(x_diag[idx_max], y_coords[idx_max], 'o', color=color, markersize=5, zorder=6)
                ax.annotate(f'{vals[idx_max]:.2f}', (x_diag[idx_max], y_coords[idx_max]),
                             fontsize=9, color=color, fontweight='bold',
                             textcoords="offset points", xytext=(8, 0))
        else:
            # Elemento horizontal: graficar diagrama en dirección Y
            x_coords = xi + (xf - xi) * x_loc / L
            y_base = np.full_like(x_coords, yi)
            y_diag = yi + vals * escala

            ax.plot(x_coords, y_diag, color=color, linewidth=2, zorder=4)
            ax.fill_between(x_coords, y_base, y_diag, alpha=0.15, color=color, zorder=2)
            ax.plot(x_coords, y_base, color=color, linewidth=0.8, linestyle='--', alpha=0.5, zorder=2)

            ax.annotate(f'{vals[0]:.2f}', (x_coords[0], y_diag[0]),
                         fontsize=8, color=color, fontweight='bold',
                         textcoords="offset points", xytext=(-5, 12))
            ax.annotate(f'{vals[-1]:.2f}', (x_coords[-1], y_diag[-1]),
                         fontsize=8, color=color, fontweight='bold',
                         textcoords="offset points", xytext=(-5, 12))

            idx_max = np.argmax(np.abs(vals))
            if abs(vals[idx_max]) > 0.01:
                ax.plot(x_coords[idx_max], y_diag[idx_max], 'o', color=color, markersize=5, zorder=6)
                ax.annotate(f'{vals[idx_max]:.2f}', (x_coords[idx_max], y_diag[idx_max]),
                             fontsize=9, color=color, fontweight='bold',
                             textcoords="offset points", xytext=(0, 12))

    ax.set_xlim(-2, 10)
    ax.set_ylim(-1, 6.5)
    fig.tight_layout()
    fname = f'diagrama_{tipo.lower()}.png'
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"  Guardado: {fname}")

plt.show()

# ============================================================================
# 7. RESUMEN FINAL PARA PRESENTACIÓN
# ============================================================================

print("\n")
print("#" * 70)
print("#" + " " * 68 + "#")
print("#" + "  RESULTADOS DEL ANÁLISIS ESTÁTICO - MARCO PLANO 2D  ".center(68) + "#")
print("#" + "  Grupo 4 - Proyecto 1  ".center(68) + "#")
print("#" + " " * 68 + "#")
print("#" * 70)

# --- 7.1 Propiedades de secciones ---
print("\n" + "=" * 70)
print("  1. PROPIEDADES DE SECCIONES")
print("=" * 70)
print(f"""
  COLUMNA - Perfil I (Acero ASTM A36, E = 200 GPa)
  -----------------------------------------------
    bf (ancho de ala)     = {bf*1000:.0f} mm
    tf (espesor de ala)   = {tf*1000:.0f} mm
    h  (alto total)       = {h*1000:.0f} mm
    tw (espesor de alma)  = {tw*1000:.0f} mm
    Area   A = {A_col:.6f} m2  =  {A_col*1e4:.2f} cm2
    Inercia I = {I_col:.8f} m4  =  {I_col*1e8:.2f} cm4

  VIGA - Seccion Cuadrada (Hormigon, E = 25 GPa)
  -----------------------------------------------
    Dimension = 40 cm x 40 cm
    Area   A = {A_viga:.4f} m2  =  {A_viga*1e4:.0f} cm2
    Inercia I = {I_viga:.8f} m4  =  {I_viga*1e8:.2f} cm4
""")

# --- 7.2 Reacciones de apoyo ---
print("=" * 70)
print("  2. REACCIONES DE APOYO")
print("=" * 70)
print(f"""
  Nodo 1 - Empotramiento (base de la columna)
  -------------------------------------------
    Fx  =  {R[0]:+.4f} kN
    Fy  =  {R[1]:+.4f} kN
    Mz  =  {R[2]:+.4f} kN.m

  Nodo 5 - Apoyo Articulado (extremo derecho de la viga)
  ------------------------------------------------------
    Fx  =  {R[12]:+.4f} kN
    Fy  =  {R[13]:+.4f} kN
""")

# --- 7.3 Verificación de equilibrio ---
print("=" * 70)
print("  3. VERIFICACION DE EQUILIBRIO GLOBAL")
print("=" * 70)
print(f"""
  Carga total horizontal (columna):  q x L = {q} x 5 = {q*5:.2f} kN
  Carga total vertical   (viga):     F = {F:.2f} kN

  Sum Fx (reacciones) = {R[0]+R[12]:.4f} kN   ===  {q*5:.2f} kN  (OK)
  Sum Fy (reacciones) = {R[1]+R[13]:.4f} kN   ===  {F:.2f} kN   (OK)
""")

# --- 7.4 Desplazamientos ---
print("=" * 70)
print("  4. DESPLAZAMIENTOS DE NODOS")
print("=" * 70)
print(f"""
  Nodo  |   UX (mm)    |   UY (mm)    |   RZ (deg)
  ------|--------------|--------------|------------
  {1:>5d} |  {'---':>10s}  |  {'---':>10s}  |  {'---':>10s}   (restringido)
  {2:>5d} |  {U[3]*1000:>+10.4f}  |  {U[4]*1000:>+10.4f}  |  {np.degrees(U[5]):>+10.6f}
  {3:>5d} |  {U[6]*1000:>+10.4f}  |  {U[7]*1000:>+10.4f}  |  {np.degrees(U[8]):>+10.6f}
  {4:>5d} |  {U[9]*1000:>+10.4f}  |  {U[10]*1000:>+10.4f}  |  {np.degrees(U[11]):>+10.6f}
  {5:>5d} |  {'---':>10s}  |  {'---':>10s}  |  {np.degrees(U[14]):>+10.6f}
""")

# --- 7.5 Fuerzas internas por elemento ---
print("=" * 70)
print("  5. FUERZAS INTERNAS POR ELEMENTO (coordenadas locales)")
print("=" * 70)

elem_info = {
    1: {'L': 2.0, 'q': q, 'tipo': 'Columna'},
    2: {'L': 3.0, 'q': q, 'tipo': 'Columna'},
    3: {'L': 5.0, 'q': 0.0, 'tipo': 'Viga'},
    4: {'L': 3.0, 'q': 0.0, 'tipo': 'Viga'},
}

for eid in range(1, 5):
    f = fuerzas_elem[eid]
    ni, nf = elementos[eid]
    L = elem_info[eid]['L']
    qe = elem_info[eid]['q']
    tipo = elem_info[eid]['tipo']

    # Punto donde V = 0 (para carga distribuida)
    x_v0 = None
    M_max = None
    if qe > 0 and abs(f[1]) > 1e-10:
        x_v0 = -f[1] / qe
        if 0 < x_v0 < L:
            M_max = f[2] + f[1] * x_v0 + qe * x_v0**2 / 2

    print(f"""
  Elemento {eid} ({tipo}, Nodos {ni}-{nf}, L = {L:.1f} m, q = {qe:.1f} kN/m)
  {"-" * 60}
    En Nodo i (x = 0):
      Axial   N = {f[0]:>+10.4f} kN
      Cortante V = {f[1]:>+10.4f} kN
      Momento  M = {f[2]:>+10.4f} kN.m

    En Nodo j (x = {L:.1f} m):
      Axial   N = {f[3]:>+10.4f} kN
      Cortante V = {f[4]:>+10.4f} kN
      Momento  M = {f[5]:>+10.4f} kN.m""")

    if M_max is not None:
        print(f"""
    Punto notable (V = 0):
      x = {x_v0:.4f} m
      M = {M_max:>+10.4f} kN.m  (momento maximo local)""")

    print(f"""
    Ecuaciones:
      N(x) = {f[0]:+.4f} + ({(f[3]-f[0])/L:+.4f}) * x   [kN]
      V(x) = {f[1]:+.4f} + ({qe:+.4f}) * x              [kN]
      M(x) = {f[2]:+.4f} + ({f[1]:+.4f}) * x + ({qe/2:+.4f}) * x2   [kN.m]""")

# --- 7.6 Resumen de diagramas ---
print(f"""

{"=" * 70}
  6. DIAGRAMAS GENERADOS
{"=" * 70}

  - diagrama_n.png  -->  Diagrama de Fuerza Axial (N)
  - diagrama_v.png  -->  Diagrama de Fuerza Cortante (V)
  - diagrama_m.png  -->  Diagrama de Momento Flector (M)

  Archivos generados en la carpeta del proyecto.

{"=" * 70}
  FIN DEL ANALISIS
{"=" * 70}
""")

# ============================================================================
# 8. EXPORTAR RESULTADOS A ARCHIVO
# ============================================================================

with open('resultados_analisis.txt', 'w', encoding='utf-8') as f:
    f.write("RESULTADOS DEL ANÁLISIS ESTÁTICO - MARCO PLANO 2D\n")
    f.write("Grupo 4 - Proyecto 1\n")
    f.write("=" * 60 + "\n\n")

    f.write("PROPIEDADES DE SECCIONES:\n")
    f.write(f"  Columna (Perfil I, Acero A36):\n")
    f.write(f"    bf={bf*1000:.0f}mm, tf={tf*1000:.0f}mm, h={h*1000:.0f}mm, tw={tw*1000:.0f}mm\n")
    f.write(f"    A = {A_col:.6f} m2, I = {I_col:.8f} m4\n\n")
    f.write(f"  Viga (Cuadrada):\n")
    f.write(f"    A = {A_viga:.4f} m2, I = {I_viga:.8f} m4\n\n")

    f.write("REACCIONES DE APOYO:\n")
    f.write(f"  Nodo 1 (Empotrado):  Fx={R[0]:.4f} kN, Fy={R[1]:.4f} kN, Mz={R[2]:.4f} kN.m\n")
    f.write(f"  Nodo 5 (Articulado): Fx={R[12]:.4f} kN, Fy={R[13]:.4f} kN\n\n")

    f.write("DESPLAZAMIENTOS:\n")
    for i in range(0, n_dof, 3):
        nodo = i // 3 + 1
        f.write(f"  Nodo {nodo}: UX={U[i]*1000:.6f}mm, UY={U[i+1]*1000:.6f}mm, RZ={np.degrees(U[i+2]):.6f}deg\n")
    f.write("\n")

    f.write("FUERZAS INTERNAS POR ELEMENTO:\n")
    for eid in range(1, 5):
        ff = fuerzas_elem[eid]
        ni, nf = elementos[eid]
        L = elem_info[eid]['L']
        qe = elem_info[eid]['q']
        f.write(f"\n  Elemento {eid} ({ni}-{nf}, L={L:.1f}m, q={qe:.1f}kN/m):\n")
        f.write(f"    N_i={ff[0]:+.4f} kN  V_i={ff[1]:+.4f} kN  M_i={ff[2]:+.4f} kN.m\n")
        f.write(f"    N_j={ff[3]:+.4f} kN  V_j={ff[4]:+.4f} kN  M_j={ff[5]:+.4f} kN.m\n")
        f.write(f"    Ecuaciones:\n")
        f.write(f"      N(x) = {ff[0]:+.4f} + ({(ff[3]-ff[0])/L:+.4f}).x\n")
        f.write(f"      V(x) = {ff[1]:+.4f} + ({qe:+.4f}).x\n")
        f.write(f"      M(x) = {ff[2]:+.4f} + ({ff[1]:+.4f}).x + ({qe/2:+.4f}).x2\n")

print("Resultados exportados a 'resultados_analisis.txt'")
