#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Section 4: The Strongly Polynomial Time Algorithm (Orlin 1993)
 *
 * Βασισμένο στο: "A Faster Strongly Polynomial Minimum Cost Flow Algorithm"
 * του James B. Orlin, Operations Research, Vol. 41, No. 2 (1993)
 *
 * Επέκταση του RHS-scaling με contractions για strongly polynomial χρόνο.
 * Διαφορές από το Section 3:
 *   1. Contraction ακμών με flow >= 3nΔ
 *   2. α = 3/4 κατώφλι για επιλογή source/sink
 *   3. Uncontraction και max-flow στο τέλος
 *
 * Μορφή αρχείου εισόδου:
 *   Γραμμή 1: n m (αριθμός κόμβων, αριθμός ακμών)
 *   Επόμενες n γραμμές: supply κάθε κόμβου
 *   Επόμενες m γραμμές: from to cost
 */

#define INF INT_MAX
#define ALPHA_NUM 3
#define ALPHA_DEN 4

typedef struct {
  int supply;           /* b(i) τρέχουσα τιμή */
  int original_supply;  /* b(i) αρχική τιμή */
  int excess;           /* e(i) πλεόνασμα */
  int potential;        /* π(i) δυναμικό */
  int distance;         /* d(i) απόσταση shortest path */
  int parent_arc;       /* γονέας στο δέντρο shortest path */
  bool active;          /* ενεργός (όχι contracted) */
} Node;

typedef struct {
  int from, to, cost, flow;
  int original_from, original_to, original_cost;
} Arc;

typedef struct {
  int node_k, node_l;
  int *potentials_snapshot;
  int n;
} ContractionRecord;

typedef struct {
  int n, m;
  Node *nodes;
  Arc *arcs;
  int delta;
  int n_active;
  ContractionRecord *history;
  int history_count;
  int history_capacity;
} Graph;

static void graph_free(Graph *g) {
  if (!g)
    return;
  for (int i = 0; i < g->history_count; i++) {
    free(g->history[i].potentials_snapshot);
  }
  free(g->nodes);
  free(g->arcs);
  free(g->history);
  free(g);
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
  g->nodes = calloc(n, sizeof(Node));
  g->arcs = calloc(m, sizeof(Arc));
  g->n_active = n;
  g->history_capacity = n;
  g->history = malloc(g->history_capacity * sizeof(ContractionRecord));
  g->history_count = 0;

  for (int i = 0; i < n; i++) {
    g->nodes[i].active = true;
    if (fscanf(fp, "%d", &g->nodes[i].supply) != 1) {
      fprintf(stderr, "Σφαλμα: Αδυναμια αναγνωσης supply κομβου %d\n", i);
      graph_free(g);
      return NULL;
    }
    g->nodes[i].original_supply = g->nodes[i].supply;
  }

  for (int i = 0; i < m; i++) {
    Arc *a = &g->arcs[i];
    if (fscanf(fp, "%d %d %d", &a->from, &a->to, &a->cost) != 3) {
      fprintf(stderr, "Σφαλμα: Αδυναμια αναγνωσης ακμης %d\n", i);
      graph_free(g);
      return NULL;
    }
    a->original_from = a->from;
    a->original_to = a->to;
    a->original_cost = a->cost;
    a->flow = 0;
  }

  return g;
}

/*
 * Αρχικοποίηση: x := 0, π := 0, e := b, Δ := 2^⌈log max{e(i)}⌉
 */
static void graph_init(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    Node *node = &g->nodes[i];
    node->excess = node->supply;
    node->potential = 0;
    node->active = true;
  }
  g->n_active = g->n;

  for (int i = 0; i < g->history_count; i++) {
    free(g->history[i].potentials_snapshot);
  }
  g->history_count = 0;

  for (int i = 0; i < g->m; i++) {
    Arc *a = &g->arcs[i];
    a->flow = 0;
    a->from = a->original_from;
    a->to = a->original_to;
    a->cost = a->original_cost;
  }

  /* Δ = 2^⌈log max{e(i)}⌉ */
  g->delta = 0;
  for (int i = 0; i < g->n; i++) {
    if (g->nodes[i].excess > g->delta)
      g->delta = g->nodes[i].excess;
  }
  if (g->delta > 0) {
    int power = 1;
    while (power < g->delta)
      power *= 2;
    g->delta = power;
  }
}

/*
 * Contraction: συγχώνευση κόμβων k και l
 * b(k) += b(l), e(k) += e(l), αντικατάσταση l με k σε όλες τις ακμές
 */
static void graph_contract(Graph *g, int k, int l) {
  if (g->history_count >= g->history_capacity) {
    g->history_capacity *= 2;
    g->history = realloc(g->history, g->history_capacity * sizeof(ContractionRecord));
  }

  ContractionRecord *rec = &g->history[g->history_count];
  rec->node_k = k;
  rec->node_l = l;
  rec->n = g->n;
  rec->potentials_snapshot = malloc(g->n * sizeof(int));
  for (int i = 0; i < g->n; i++)
    rec->potentials_snapshot[i] = g->nodes[i].potential;
  g->history_count++;

  /* Τα reduced costs γίνονται νέα κόστη ακμών */
  for (int i = 0; i < g->m; i++) {
    int from = g->arcs[i].from, to = g->arcs[i].to;
    if (from != to)
      g->arcs[i].cost -= g->nodes[from].potential - g->nodes[to].potential;
  }

  for (int i = 0; i < g->n; i++)
    g->nodes[i].potential = 0;

  g->nodes[k].supply += g->nodes[l].supply;
  g->nodes[k].excess += g->nodes[l].excess;

  for (int i = 0; i < g->m; i++) {
    if (g->arcs[i].from == l) g->arcs[i].from = k;
    if (g->arcs[i].to == l)   g->arcs[i].to = k;
  }

  g->nodes[l].active = false;
  g->n_active--;
}

/*
 * Contraction ακμών με x_ij >= 3nΔ (strongly feasible)
 */
static void graph_contract_strongly_feasible(Graph *g) {
  bool contracted;
  do {
    contracted = false;
    for (int i = 0; i < g->m; i++) {
      int from = g->arcs[i].from, to = g->arcs[i].to;
      if (from == to) continue;
      if (!g->nodes[from].active || !g->nodes[to].active) continue;

      if (g->arcs[i].flow >= 3LL * g->n * g->delta) {
        graph_contract(g, from, to);
        contracted = true;
        break;
      }
    }
  } while (contracted);
}

/*
 * S(Δ) = {i : e(i) >= αΔ} με α = 3/4
 */
static int graph_find_source(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    Node *node = &g->nodes[i];
    if (node->active && (long long)ALPHA_DEN * node->excess >= (long long)ALPHA_NUM * g->delta)
      return i;
  }
  return -1;
}

/*
 * T(Δ) = {i : e(i) <= -αΔ} με α = 3/4
 */
static int graph_find_sink(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    Node *node = &g->nodes[i];
    if (node->active && (long long)ALPHA_DEN * node->excess <= -(long long)ALPHA_NUM * g->delta)
      return i;
  }
  return -1;
}

/*
 * Αλγόριθμος Dijkstra στο residual network G(x)
 */
static void graph_compute_shortest_paths(Graph *g, int source) {
  bool *visited = calloc(g->n, sizeof(bool));

  for (int i = 0; i < g->n; i++) {
    g->nodes[i].distance = INF;
    g->nodes[i].parent_arc = -1;
  }
  g->nodes[source].distance = 0;

  for (int round = 0; round < g->n; round++) {
    int choice = -1, min_distance = INF;
    for (int i = 0; i < g->n; i++) {
      if (g->nodes[i].active && !visited[i] && g->nodes[i].distance < min_distance) {
        min_distance = g->nodes[i].distance;
        choice = i;
      }
    }
    if (choice == -1) break;
    visited[choice] = true;

    int distance_choice = g->nodes[choice].distance;
    for (int i = 0; i < g->m; i++) {
      int from = g->arcs[i].from, to = g->arcs[i].to;
      if (from == to) continue;

      int rc = g->arcs[i].cost - g->nodes[from].potential + g->nodes[to].potential;

      /* Forward arc */
      if (from == choice && g->nodes[to].active && !visited[to]) {
        int new_distance = distance_choice + rc;
        if (new_distance < g->nodes[to].distance) {
          g->nodes[to].distance = new_distance;
          g->nodes[to].parent_arc = i;
        }
      }
      /* Reverse arc (μόνο αν flow > 0) */
      if (to == choice && g->arcs[i].flow > 0 && g->nodes[from].active && !visited[from]) {
        int new_distance = distance_choice - rc;
        if (new_distance < g->nodes[from].distance) {
          g->nodes[from].distance = new_distance;
          g->nodes[from].parent_arc = -(i + 1);
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
    Node *node = &g->nodes[i];
    if (node->active && node->distance != INF)
      node->potential -= node->distance;
  }
}

/*
 * Αύξηση (augment) Δ μονάδων κατά μήκος του shortest path από source σε sink
 * Ενημέρωση πλεονασμάτων: e(source) -= Δ, e(sink) += Δ
 */
static void graph_augment_flow(Graph *g, int source, int sink) {
  int current = sink;
  while (current != source) {
    int parent = g->nodes[current].parent_arc;
    if (parent >= 0) {
      g->arcs[parent].flow += g->delta;
      current = g->arcs[parent].from;
    } else {
      int reverse = -(parent + 1);
      g->arcs[reverse].flow -= g->delta;
      current = g->arcs[reverse].to;
    }
  }
  g->nodes[source].excess -= g->delta;
  g->nodes[sink].excess += g->delta;
}

/*
 * Uncontraction: αναίρεση contractions με αντίστροφη σειρά
 */
static void graph_uncontract_all(Graph *g) {
  for (int h = g->history_count - 1; h >= 0; h--) {
    ContractionRecord *rec = &g->history[h];
    int k = rec->node_k, l = rec->node_l;

    g->nodes[l].active = true;
    g->n_active++;
    g->nodes[l].potential = g->nodes[k].potential;

    for (int i = 0; i < g->n; i++)
      g->nodes[i].potential += rec->potentials_snapshot[i];
  }

  for (int i = 0; i < g->m; i++) {
    g->arcs[i].from = g->arcs[i].original_from;
    g->arcs[i].to = g->arcs[i].original_to;
  }
}

/*
 * Section 4.5: Max-flow στο zero-reduced-cost subnetwork G°
 *
 * Αλγόριθμος (από το paper):
 *   1. G° = {(i,j) ∈ A : c̄ij = 0} όπου c̄ij = cij - π(i) + π(j)
 *   2. Προσθήκη s* με ακμές προς κόμβους supply (χωρητικότητα b(i))
 *   3. Προσθήκη t* με ακμές από κόμβους demand (χωρητικότητα -b(i))
 *   4. Max-flow από s* προς t*
 */

/* Τύποι ακμών στο augmenting path */
enum { EDGE_SOURCE, EDGE_FORWARD, EDGE_REVERSE, EDGE_SINK };

/*
 * BFS: Εύρεση augmenting path στο G° ∪ {s*, t*}
 */
static bool maxflow_bfs(Graph *g, int *residual_src, int *residual_sink,
                        int *parent, int *edge_type, int *edge_arc,
                        bool *visited, int *queue) {
  const int SUPER_SRC = g->n;
  const int SUPER_SINK = g->n + 1;
  const int total_nodes = g->n + 2;

  /* Αρχικοποίηση */
  for (int i = 0; i < total_nodes; i++) {
    visited[i] = false;
    parent[i] = -1;
  }

  /* Ξεκινάμε BFS από τον SUPER_SRC */
  int front = 0, back = 0;
  queue[back++] = SUPER_SRC;
  visited[SUPER_SRC] = true;

  while (front < back) {
    int current = queue[front++];

    /* Από τον SUPER_SRC, προσθέτουμε στην ουρά κόμβους με εναπομείναν supply */
    if (current == SUPER_SRC) {
      for (int node = 0; node < g->n; node++) {
        if (residual_src[node] > 0) {
          visited[node] = true;
          parent[node] = SUPER_SRC;
          edge_type[node] = EDGE_SOURCE;
          queue[back++] = node;
        }
      }
      continue;
    }

    /* Αν ο τρέχων κόμβος έχει εναπομείναν demand, βρήκαμε μονοπάτι */
    if (residual_sink[current] > 0) {
      parent[SUPER_SINK] = current;
      edge_type[SUPER_SINK] = EDGE_SINK;
      return true;
    }

    /* Ακολουθούμε ακμές μηδενικού reduced cost στο G° */
    for (int arc = 0; arc < g->m; arc++) {
      int from = g->arcs[arc].original_from;
      int to = g->arcs[arc].original_to;
      int rc = g->arcs[arc].original_cost - g->nodes[from].potential + g->nodes[to].potential;

      if (rc != 0)
        continue;

      /* Forward arc: current → to */
      if (from == current && !visited[to]) {
        visited[to] = true;
        parent[to] = current;
        edge_type[to] = EDGE_FORWARD;
        edge_arc[to] = arc;
        queue[back++] = to;
      }

      /* Reverse arc: current → from (μόνο αν flow > 0) */
      if (to == current && g->arcs[arc].flow > 0 && !visited[from]) {
        visited[from] = true;
        parent[from] = current;
        edge_type[from] = EDGE_REVERSE;
        edge_arc[from] = arc;
        queue[back++] = from;
      }
    }
  }

  return false;
}

/*
 * Υπολογισμός βέλτιστης ροής με Ford-Fulkerson (bottleneck augmentation)
 */
static void graph_compute_optimal_flow(Graph *g) {
  /* Μηδενισμός ροών */
  for (int i = 0; i < g->m; i++)
    g->arcs[i].flow = 0;

  /* Δέσμευση βοηθητικών πινάκων */
  int *residual_src = calloc(g->n, sizeof(int));   /* χωρητικότητα s* → v */
  int *residual_sink = calloc(g->n, sizeof(int));  /* χωρητικότητα v → t* */
  int *parent = malloc((g->n + 2) * sizeof(int));
  int *edge_type = malloc((g->n + 2) * sizeof(int));
  int *edge_arc = malloc((g->n + 2) * sizeof(int));
  bool *visited = malloc((g->n + 2) * sizeof(bool));
  int *queue = malloc((g->n + 2) * sizeof(int));

  /* Αρχικοποίηση χωρητικοτήτων από supplies/demands */
  for (int i = 0; i < g->n; i++) {
    int supply = g->nodes[i].original_supply;
    if (supply > 0)      residual_src[i] = supply;
    else if (supply < 0) residual_sink[i] = -supply;
  }

  const int SUPER_SRC = g->n, SUPER_SINK = g->n + 1;

  /* Ford-Fulkerson με augmentation κατά bottleneck */
  while (maxflow_bfs(g, residual_src, residual_sink,
                     parent, edge_type, edge_arc, visited, queue)) {

    /* Βήμα 1: Εύρεση bottleneck — ελάχιστη χωρητικότητα κατά μήκος μονοπατιού */
    int bottleneck = INF;
    int step = SUPER_SINK;
    while (step != SUPER_SRC) {
      int capacity = INF;

      if (edge_type[step] == EDGE_SOURCE)
        capacity = residual_src[step];
      else if (edge_type[step] == EDGE_SINK)
        capacity = residual_sink[parent[step]];
      else if (edge_type[step] == EDGE_REVERSE)
        capacity = g->arcs[edge_arc[step]].flow;
      /* EDGE_FORWARD: απεριόριστη χωρητικότητα στο G° */

      if (capacity < bottleneck)
        bottleneck = capacity;

      step = parent[step];
    }

    if (bottleneck <= 0 || bottleneck == INF)
      break;

    /* Βήμα 2: Αύξηση ροής κατά bottleneck κατά μήκος μονοπατιού */
    step = SUPER_SINK;
    while (step != SUPER_SRC) {
      if (edge_type[step] == EDGE_SOURCE)
        residual_src[step] -= bottleneck;
      else if (edge_type[step] == EDGE_SINK)
        residual_sink[parent[step]] -= bottleneck;
      else if (edge_type[step] == EDGE_FORWARD)
        g->arcs[edge_arc[step]].flow += bottleneck;
      else if (edge_type[step] == EDGE_REVERSE)
        g->arcs[edge_arc[step]].flow -= bottleneck;

      step = parent[step];
    }
  }

  free(residual_src);
  free(residual_sink);
  free(parent);
  free(edge_type);
  free(edge_arc);
  free(visited);
  free(queue);
}

static bool has_imbalanced_node(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    if (g->nodes[i].active && g->nodes[i].excess != 0)
      return true;
  }
  return false;
}

/*
 * Βοηθητικές συναρτήσεις complete regeneration (Figure 5)
 * Έλεγχος αν xij = 0 για όλες τις ακμές στο contracted γράφο
 */
static bool all_flows_zero(Graph *g) {
  for (int i = 0; i < g->m; i++) {
    if (g->arcs[i].from != g->arcs[i].to && g->arcs[i].flow != 0)
      return false;
  }
  return true;
}

/* Έλεγχος αν e(i) < Δ για όλους τους ενεργούς κόμβους */
static bool all_excesses_below_delta(Graph *g) {
  for (int i = 0; i < g->n; i++) {
    if (g->nodes[i].active && g->nodes[i].excess >= g->delta)
      return false;
  }
  return true;
}

/* Υπολογισμός max{e(i)} για ενεργούς κόμβους με θετικό πλεόνασμα */
static int max_positive_excess(Graph *g) {
  int max_excess = 0;
  for (int i = 0; i < g->n; i++) {
    if (g->nodes[i].active && g->nodes[i].excess > max_excess)
      max_excess = g->nodes[i].excess;
  }
  return max_excess;
}

/*
 * Κύριος αλγόριθμος Section 4
 */
static bool graph_run_main_loop(Graph *g) {
  while (has_imbalanced_node(g)) {
    /*
     * Complete regeneration (Figure 5):
     * if xij = 0 for all (i,j) ∈ Â and e(i) < Δ for all i ∈ N̂ then
     *     Δ := max{e(i): i ∈ N̂}
     */
    if (all_flows_zero(g) && all_excesses_below_delta(g)) {
      int new_delta = max_positive_excess(g);
      if (new_delta == 0) break;  /* Όλοι οι κόμβοι ισορροπημένοι */
      int power = 1;
      while (power < new_delta) power *= 2;
      g->delta = power;
    }

    graph_contract_strongly_feasible(g);

    while (true) {
      int source = graph_find_source(g);
      int sink = graph_find_sink(g);
      if (source == -1 || sink == -1) break;

      graph_compute_shortest_paths(g, source);
      if (g->nodes[sink].distance == INF) return false;

      graph_update_potentials(g);
      graph_augment_flow(g, source, sink);
    }

    g->delta /= 2;

    if (g->delta == 0 && has_imbalanced_node(g))
      break;
  }

  graph_uncontract_all(g);
  graph_compute_optimal_flow(g);

  return true;
}

/*
 * Συνολικό κόστος: Σ c_ij * x_ij
 */
static long long graph_get_total_cost(Graph *g) {
  long long total = 0;
  for (int i = 0; i < g->m; i++) {
    total += (long long)g->arcs[i].flow * g->arcs[i].original_cost;
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
