# Análisis Estático de Marco Plano 2D - Hiperestático

**Grupo 4 - Proyecto 1**

## Descripción del Problema

Análisis estático de un marco plano 2D hiperestático compuesto por una columna vertical empotrada en la base y una viga horizontal articulada en el extremo derecho.

### Geometría

```
Nodo 3 (0,5)  ← Tope columna (libre)
    │
    │  Columna (Perfil I, Acero A36)
    │  q = 17 kN/m →
    │
Nodo 2 (0,2)  ← Unión rígida columna-viga
    │─────────────────────────●─────────── Nodo 5 (8,2)
    │    Viga (40x40 cm)     F=20kN↓      Apoyo articulado
    │                         Nodo 4 (5,2)
    │
Nodo 1 (0,0)  ← Empotramiento
```

### Dimensiones
- **Columna:** 5 m de altura (2 m + 3 m)
- **Viga:** 8 m de longitud (5 m + 3 m)

### Cargas
- **Columna:** Carga distribuida horizontal q = 17 kN/m (dirección +X)
- **Viga:** Carga puntual vertical F = 20 kN (dirección -Y, en X = 5 m)

### Condiciones de Borde
- **Nodo 1:** Empotramiento (restringe UX, UY, RZ)
- **Nodo 5:** Articulación (restringe UX, UY, permite RZ)

---

## Propiedades de Secciones

### Columna - Perfil I (Acero ASTM A36, E = 200 GPa)

| Dimensión | Valor |
|-----------|-------|
| bf (ancho de ala) | 200 mm |
| tf (espesor de ala) | 12 mm |
| h (alto total) | 300 mm |
| tw (espesor de alma) | 8 mm |

- **Área (A):** 57.60 cm²
- **Inercia (I):** 52,546.40 cm⁴

### Viga - Sección Cuadrada (Hormigón, E = 25 GPa)

| Dimensión | Valor |
|-----------|-------|
| Base | 40 cm |
| Alto | 40 cm |

- **Área (A):** 1,600 cm²
- **Inercia (I):** 213,333.33 cm⁴

---

## Método de Solución

El script utiliza dos métodos independientes para verificación cruzada:

### 1. Método Matricial Directo
- Ensambla la matriz de rigidez global del sistema
- Aplica las cargas equivalentes nodales (incluyendo momentos equivalentes)
- Resuelve el sistema de ecuaciones K·U = F
- Calcula reacciones y fuerzas internas

### 2. OpenSeesPy
- Análisis computacional con el software OpenSees
- Permite verificar los resultados del método matricial

---

## Archivos del Repositorio

| Archivo | Descripción |
|---------|-------------|
| `marco_2d_asentamiento.py` | Script principal con ambos métodos de análisis |
| `diagrama_n.png` | Diagrama de fuerza axial |
| `diagrama_v.png` | Diagrama de fuerza cortante |
| `diagrama_m.png` | Diagrama de momento flector |
| `resultados_analisis.txt` | Resumen de resultados numéricos |

---

## Requisitos

```bash
pip install numpy matplotlib
pip install openseespy  # Opcional, para verificación
```

---

## Ejecución

```bash
python marco_2d_asentamiento.py
```

---

## Resultados Principales

### Reacciones de Apoyo

| Nodo | Tipo | Fx (kN) | Fy (kN) | Mz (kN·m) |
|------|------|---------|---------|-----------|
| 1 | Empotrado | *Ver salida* | *Ver salida* | *Ver salida* |
| 5 | Articulado | *Ver salida* | *Ver salida* | - |

### Verificación de Equilibrio

```
ΣFx = Rx1 + Rx5 = 85 kN (carga total: 85 kN)
ΣFy = Ry1 + Ry5 = 20 kN (carga total: 20 kN)
```

---

## Ecuaciones de Fuerzas Internas

Para cada elemento, las fuerzas internas se expresan en función de la coordenada local x:

### Columna (Elementos 1 y 2) - Con carga distribuida q = 17 kN/m

```
N(x) = Ni + (Nj - Ni)·x/L          [kN]
V(x) = Vi + q·x                     [kN]
M(x) = Mi + Vi·x + q·x²/2          [kN·m]
```

### Viga (Elementos 3 y 4) - Sin carga distribuida

```
N(x) = Ni + (Nj - Ni)·x/L          [kN]
V(x) = Vi                           [kN]  (constante)
M(x) = Mi + Vi·x                    [kN·m]  (lineal)
```

---

## Diagramas de Fuerzas Internas

Los diagramas se generan superpuestos a la geometría real de la estructura:

- **Fuerza Axial (N):** Mostrada en rojo, perpendicular al elemento
- **Fuerza Cortante (V):** Mostrada en azul, perpendicular al elemento
- **Momento Flector (M):** Mostrado en verde, perpendicular al elemento

Los valores positivos se grafican en la dirección positiva del eje local del elemento.

---

## Créditos

**Grupo 4 - Proyecto 1**
Ingeniero Estructural
