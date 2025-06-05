# GbSMTWTP: Operadores de caja gris para SMTWTP y algoritmo DRILS

## Características
- **DRILS**: Clase que implementa el algoritmo DRILS y el operador de cruce *Partition Crossover*
- **Neighborhood**: Implementación de vecindario para búsqueda local
- **local_search**: Operador de búsqueda local para SMTWTP
- **SMWTP**: Implementación de la función objetivo y delta-evaluación.

## Instalación

### Requisitos del sistema

- Linux (soporte principal y probado)
- C++17
- [CMake](https://cmake.org/) >= 3.10, `make`

### Comandos para instalación

```bash
git clone https://github.com/lvargasgarcia/GbSMTWTP.git
cd GbSMTWTP
chmod 700 install.sh
sudo ./install.sh

```
Esto moverá los archivos *.h a /usr/local/include/
---

### Funciones principales

#### DRILS
```cpp
std::vector<std::tuple<long, long, double>> basic_DRILS(SMWTP& instance, int neighborhood_size, std::mt19937& device, long time_to_run, double perturbation_factor, int tolerance)
```
- **Descripción:**  
  Ejecuta el algoritmo DRILS sobre una instancia de SMWTP (Single Machine Weighted Tardiness Problem) durante un tiempo dado, aplicando búsqueda local y perturbaciones controladas.
- **Parámetros:**
  - `instance`: Referencia a un objeto de tipo `SMWTP` que representa el problema.
  - `neighborhood_size`: Tamaño del vecindario para la búsqueda local.
  - `device`: Generador de números aleatorios (`std::mt19937`).
  - `time_to_run`: Tiempo máximo de ejecución (en segundos).
  - `perturbation_factor`: Swaps en cada iteración (0.0 a 1.0) (se efectúan perturbation_factor*N) swaps.
  - `tolerance`: Número de iteraciones sin mejora antes de aceptar una nueva solución.
- **Retorno:**  
  Vector de tuplas `(tiempo_ms, mejor_valor, porcentaje_px)` con el historial de mejoras y porcentaje de éxito en el operador de cruce.

#### Partition crossover 
```cpp
permutation_t partition_crossover(const permutation_t& sigma_1, const permutation_t& sigma_2, SMWTP& instance)
```

**Descripción:**  
Realiza el operador de cruce *Partition Crossover* (PX) entre dos permutaciones, generando un hijo que combina bloques de ambos padres, optimizando localmente según la función objetivo de la instancia.

**Parámetros:**
- `sigma_1`: Primera permutación padre.
- `sigma_2`: Segunda permutación padre.
- `instance`: Referencia a la instancia de `SMWTP` para evaluar los bloques y calcular mejoras.

**Retorno:**  
Una permutación (`permutation_t`) que representa el hijo resultante del cruce PX entre los padres.

---
```cpp
std::pair<permutation_t, long> local_search(Neighborhood neighborhood, SMWTP& instance, const permutation_t& pi, std::mt19937& device)
```

**Descripción:**  
Aplica búsqueda local iterativa sobre una permutación inicial, explorando el vecindario definido y aceptando movimientos que mejoran la función objetivo.

**Parámetros:**
- `neighborhood`: Objeto `Neighborhood` que define los movimientos permitidos en el vecindario.
- `instance`: Referencia a la instancia de `SMWTP` para evaluar soluciones.
- `pi`: Permutación inicial sobre la que se realiza la búsqueda local.
- `device`: Generador de números aleatorios (`std::mt19937`) para seleccionar movimientos.

**Retorno:**  
Un par (`std::pair<permutation_t, long>`) donde:
- El primer elemento es la permutación óptima local encontrada.
- El segundo elemento es el valor de la función objetivo asociado a esa permutación.

#### Clase Neighborhood

## Clase `Neighborhood`

La clase [`Neighborhood`](include/Neighborhood.hpp) define y gestiona el vecindario de movimientos para la búsqueda local en SMTWTP.

### Constructores

```cpp
Neighborhood(int n, int k);
```
- **Parámetros:**
  - `n`: Tamaño de la permutación (número de trabajos).
  - `k`: Tamaño de los bloques de permutación considerados en el vecindario (por ejemplo, para bloques de swaps de tamaño `k`).

### Atributos principales

- `perms_set_t perms_set`: Conjunto de permutaciones posibles para bloques de tamaño `k`.
- `int n`: Tamaño de la instancia asociada.
- `int k`: Elementos consecutivos movidos por los movimientos considerados.
- `long size`: Número total de movimientos posibles en el vecindario.
- `int categories`: Número de categorías de bloques (n - k + 1).
- `long perms_count`: Número de permutaciones posibles para un bloque de tamaño `k` (factorial(k) - 1).
- `std::vector<std::pair<int,int>> improving_moves_set`: Lista de movimientos que mejoran la solución actual.
- `std::vector<int> improving_moves_indexes`: Índices para acceso rápido a los movimientos en `improving_moves_set`.
- `std::vector<long> scores`: Mejora (delta) asociada a cada movimiento.
- `std::vector<long> common_tardiness_set`: Información auxiliar para la evaluación de movimientos.

### Métodos principales

#### `insert_move(const int position, const int index)`

Agrega un movimiento a la lista de movimientos que mejoran la solución, si aún no está presente.

- **Parámetros:**
  - `position`: Posición del bloque en la permutación.
  - `index`: Índice del movimiento en el conjunto de movimientos posibles.

#### `remove_move(const int index)`

Elimina un movimiento de la lista de movimientos que mejoran la solución.

- **Parámetros:**
  - `index`: Índice del movimiento a eliminar.

#### `std::pair<int,int> select_improving_move(std::mt19937& gen)`

Selecciona aleatoriamente un movimiento de entre los que mejoran la solución.

- **Parámetros:**
  - `gen`: Generador de números aleatorios.
- **Retorno:**  
  Par `(position, index)` que identifica el movimiento seleccionado.

