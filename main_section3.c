#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Section 3: The Edmonds-Karp Scaling Technique (Orlin 1993)
 *
 * Βασισμένο στο: "A Faster Strongly Polynomial Minimum Cost Flow Algorithm"
 * του James B. Orlin, Operations Research, Vol. 41, No. 2 (1993)
 *
 * Μορφή αρχείου εισόδου:
 *   Γραμμή 1: n m (αριθμός κόμβων, αριθμός ακμών)
 *   Επόμενες n γραμμές: supply κάθε κόμβου
 *   Επόμενες m γραμμές: from to cost
 */

#define INF INT_MAX

typedef struct {
  int from;
  int to;
  int cost;
  int flow;
} Arc;

typedef struct {
  int n;
  int m;
  int *supplies;
  int *excess;
  int *potentials;
  int *distances;
  int *parent_arc;
  int *reduced_costs;
  Arc *arcs;
  int delta;
} Graph;

static void graph_free(Graph *g) {
  if (!g)
    return;
  free(g->supplies);
  free(g->excess);
  free(g->potentials);
  free(g->distances);
  free(g->parent_arc);
  free(g->reduced_costs);
  free(g->arcs);
  free(g);
}

static int ceil_log2(int x) {
  if (x <= 1)
    return 0;
  int log = 0;
  int val = 1;
  while (val < x) {
    val *= 2;
    log++;
  }
  return log;
}


static Graph *graph_create(FILE *fp) {
  int n, m;
  if (fscanf(fp, "%d %d", &n, &m) != 2) {
    fprintf(stderr, "Σφαλμα: Μη εγκυρη μορφη αρχειου\n");
    return NULL;
  }

  Graph *g = malloc(sizeof(Graph));
  g->n = n;
  g->m = m;
  g->supplies = calloc(n, sizeof(int));
  g->excess = calloc(n, sizeof(int));
  g->potentials = calloc(n, sizeof(int));
  g->distances = calloc(n, sizeof(int));
  g->parent_arc = malloc(n * sizeof(int));
  g->reduced_costs = malloc(m * sizeof(int));
  g->arcs = calloc(m, sizeof(Arc));

  for (int i = 0; i < n; i++) {
    if (fscanf(fp, "%d", &g->supplies[i]) != 1) {
      fprintf(stderr, "Σφαλμα: Αδυναμια αναγνωσης supply κομβου %d\n", i);
      graph_free(g);
      return NULL;
    }
  }

  for (int i = 0; i < m; i++) {
    if (fscanf(fp, "%d %d %d", &g->arcs[i].from, &g->arcs[i].to,
               &g->arcs[i].cost) != 3) {
      fprintf(stderr, "Σφαλμα: Αδυναμια αναγνωσης ακμης %d\n", i);
      graph_free(g);
      return NULL;
    }
    g->arcs[i].flow = 0;
  }

  return g;
}

/*
 * Αρχικοποίηση: x := 0, π := 0, e := b
 * U = max{|b(i)| : i ∈ N}
 * Δ = 2^⌈log U⌉
 */
static void graph_init(Graph *g) {
  int U = 0;
  for (int i = 0; i < g->n; i++) {
    int abs_supply = abs(g->supplies[i]);
    if (abs_supply > U)
      U = abs_supply;
  }

  g->delta = (U == 0) ? 1 : (1 << ceil_log2(U));

  for (int i = 0; i < g->n; i++) {
    g->excess[i] = g->supplies[i];
    g->potentials[i] = 0;
  }

  for (int i = 0; i < g->m; i++) {
    g->arcs[i].flow = 0;
  }
}

/*
 * Μειωμένο κόστος (Reduced cost): c^π_ij = c_ij - π(i) + π(j)
 */
static void graph_compute_reduced_costs(Graph *g) {
  for (int i = 0; i < g->m; i++) {
    int from = g->arcs[i].from;
    int to = g->arcs[i].to;
    g->reduced_costs[i] =
        g->arcs[i].cost - g->potentials[from] + g->potentials[to];
  }
}

/*
 * Αλγόριθμος Dijkstra στο residual network G(x)
 */
static void graph_compute_shortest_paths(Graph *g, int source) {
  bool *visited = calloc(g->n, sizeof(bool));

  for (int i = 0; i < g->n; i++) {
    g->distances[i] = INF;
    g->parent_arc[i] = -1;
  }
  g->distances[source] = 0;

  for (int round = 0; round < g->n; round++) {
    /* Βρες τον μη επισκεφθέντα κόμβο με ελάχιστη απόσταση */
    int choice = -1;
    int min_distance = INF;
    for (int i = 0; i < g->n; i++) {
      if (!visited[i] && g->distances[i] < min_distance) {
        min_distance = g->distances[i];
        choice = i;
      }
    }

    if (choice == -1)
      break;
    visited[choice] = true;

    /* Χαλάρωση (relaxation) όλων των ακμών από τον κόμβο choice */
    for (int i = 0; i < g->m; i++) {
      int from = g->arcs[i].from;
      int to = g->arcs[i].to;

      /* Forward arc: choice -> to */
      if (from == choice && !visited[to]) {
        int new_distance = g->distances[choice] + g->reduced_costs[i];
        if (new_distance < g->distances[to]) {
          g->distances[to] = new_distance;
          g->parent_arc[to] = i;
        }
      }

      /* Reverse arc: choice -> from (μόνο αν flow > 0) */
      if (to == choice && g->arcs[i].flow > 0 && !visited[from]) {
        int new_distance = g->distances[choice] - g->reduced_costs[i];
        if (new_distance < g->distances[from]) {
          g->distances[from] = new_distance;
          g->parent_arc[from] = -(i + 1);
        }
      }
    }
  }

  free(visited);
}

/*
 * Ενημέρωση potentials: π(i) := π(i) - d(i)
 */
static void graph_update_potentials(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    if (g->distances[i] != INF) {
      g->potentials[i] -= g->distances[i];
    }
  }
}

/*
 * Αύξηση (augment) Δ μονάδων κατά μήκος του shortest path από source σε sink
 * Ενημέρωση πλεονασμάτων: e(source) -= Δ, e(sink) += Δ
 */
static void graph_augment_flow(Graph *g, int source, int sink) {
  int current = sink;

  while (current != source) {
    int parent = g->parent_arc[current];

    if (parent >= 0) {
      g->arcs[parent].flow += g->delta;
      current = g->arcs[parent].from;
    } else {
      int reverse = -(parent + 1);
      g->arcs[reverse].flow -= g->delta;
      current = g->arcs[reverse].to;
    }
  }

  g->excess[source] -= g->delta;
  g->excess[sink] += g->delta;
}

/*
 * S(Δ) = {i : e(i) >= Δ} - Σύνολο κόμβων με πλεόνασμα (excess)
 */
static int graph_find_source(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    if (g->excess[i] >= g->delta)
      return i;
  }
  return -1;
}

/*
 * T(Δ) = {i : e(i) <= -Δ} - Σύνολο κόμβων με έλλειμμα (deficit)
 */
static int graph_find_sink(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    if (g->excess[i] <= -g->delta)
      return i;
  }
  return -1;
}

/*
 * Κύριος αλγόριθμος Section 3
 */
static bool graph_run_main_loop(Graph *g) {
  while (g->delta >= 1) {
    while (true) {
      int source = graph_find_source(g);
      int sink = graph_find_sink(g);

      if (source == -1 || sink == -1)
        break;

      graph_compute_reduced_costs(g);
      graph_compute_shortest_paths(g, source);

      if (g->distances[sink] == INF)
        return false;

      graph_update_potentials(g);
      graph_augment_flow(g, source, sink);
    }

    g->delta /= 2;
  }

  return true;
}

/*
 * Συνολικό κόστος: Σ c_ij * x_ij
 */
static long long graph_get_total_cost(Graph *g) {
  long long total = 0;
  for (int i = 0; i < g->m; i++) {
    total += (long long)g->arcs[i].flow * g->arcs[i].cost;
  }
  return total;
}

int main(void) {
  Graph *g = graph_create(stdin);
  if (!g)
    return 1;

  graph_init(g);

  bool success = graph_run_main_loop(g);

  if (success) {
    printf("Βελτιστο κοστος: %lld\n", graph_get_total_cost(g));
  } else {
    printf("Το προβλημα δεν εχει εφικτη λυση\n");
  }

  graph_free(g);
  return success ? 0 : 1;
}
